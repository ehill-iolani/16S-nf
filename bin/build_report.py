#!/usr/bin/env python3
"""
Aggregate per-cluster BLAST hits + consensus fastas into a single
sample x taxon abundance table, plus a minimal QC/summary HTML.

This is a functional starting point, not the final report -- extend with
per-cluster read counts (from vsearch clusters.uc), top-hit
filtering by pident/evalue, and a proper MultiQC-style layout once the
sample sheet format and reference taxonomy fields are finalized.
"""
import argparse
import os
import re
import shutil

import pandas as pd

COLS = ["seq_id", "subject_id", "pident", "length", "evalue", "bitscore", "stitle"]

# SILVA's export puts a lineage where other references put a binomial: the
# description BLAST reports for a SILVA hit is
#   "Bacteria;Proteobacteria;Gammaproteobacteria;...;Pseudomonas;<organism>"
# SILVA's curated taxonomy stops at genus -- the trailing organism name is just
# whatever the source record said ("uncultured bacterium" for most of them) --
# so a SILVA hit is called at genus level, and the lineage is written out as
# rank columns named like edna-ont-nf's --taxdump ones, which the frontend
# already rolls up by. SILVA's "domain" goes in NCBI's "superkingdom" column.
LINEAGE_COLS = ["superkingdom", "phylum", "class", "order", "family", "genus"]
SILVA_RANK_COL = {"domain": "superkingdom", "phylum": "phylum", "class": "class",
                  "order": "order", "family": "family", "genus": "genus"}
# Names SILVA uses for a rank it couldn't fill; treated as no name at all, so
# one "uncultured" family doesn't lump unrelated clusters together.
SILVA_PLACEHOLDER = re.compile(r"^(uncultured|incertae sedis|unknown\b)", re.I)


def is_silva(stitle):
    return isinstance(stitle, str) and ";" in stitle


def load_silva_taxmap(path):
    """SILVA's tax_slv_ssu file: "path;\\ttaxid\\trank\\t..." -> {path: rank}."""
    taxmap = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                taxmap[fields[0]] = fields[2]
    return taxmap


def silva_lineage(stitle, taxmap=None):
    """{rank column: name} for a SILVA lineage description. Each level's rank
    comes from SILVA's own rank map when there is one; without it (a
    SILVA-format --taxdb from somewhere else), levels are taken in the usual
    domain > phylum > ... > genus order, which is right for all but ~0.4% of
    SILVA 138.1 records (the ones that skip a rank). Placeholders are left out."""
    names = stitle.split(";")[:-1]  # the last field is the organism name
    out = {}
    for i, name in enumerate(names):
        rank = taxmap.get(";".join(names[:i + 1]) + ";") if taxmap else None
        col = SILVA_RANK_COL.get(rank) if rank else (LINEAGE_COLS[i] if i < len(LINEAGE_COLS) else None)
        name = name.strip()
        if col and name and not SILVA_PLACEHOLDER.match(name):
            out[col] = name
    return out


def lineage_label(lineage):
    """A SILVA hit's taxon label: "Genus sp.", or when the genus is unnamed, the
    lowest named rank with that rank in brackets ("Anaerolineaceae (family)").
    Mirrored by pcoa.py and the frontend's resultsFormat.js -- keep in sync."""
    if lineage.get("genus"):
        return f"{lineage['genus']} sp."
    for col in reversed(LINEAGE_COLS[:-1]):
        if lineage.get(col):
            return f"{lineage[col]} ({'domain' if col == 'superkingdom' else col})"
    return "Unclassified"


def load_hits(paths):
    frames = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t", names=COLS)
            df["source_file"] = p
            frames.append(df)
        except pd.errors.EmptyDataError:
            continue
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLS)


def species_label(stitle, taxmap=None):
    # same rule as the frontend's shortTaxonLabel: BLAST's stitle starts with
    # the binomial, so the first two tokens name the species -- except for a
    # SILVA lineage, which is called at genus level (see lineage_label)
    if is_silva(stitle):
        return lineage_label(silva_lineage(stitle, taxmap))
    toks = stitle.split() if isinstance(stitle, str) else []
    return " ".join(toks[:2])


def tied_species(hits_df, taxmap=None):
    """seq_id -> "|"-joined distinct species among the hits tied for that
    cluster's best bitscore, or "" when the best hit is unambiguous.

    A cluster's hits are often several reference records with identical
    bitscores; when they name different species the working call (the first
    of them) is really a coin flip, and the frontend reports it at genus
    level (same genus) or as ambiguous (different genera) instead."""
    out = {}
    hits = hits_df[hits_df["subject_id"] != "NO_HIT"]
    for seq_id, grp in hits.groupby("seq_id", sort=False):
        top = grp[grp["bitscore"] == grp["bitscore"].max()]
        names = list(dict.fromkeys(n for n in (species_label(s, taxmap) for s in top["stitle"]) if n))
        out[seq_id] = "|".join(names) if len(names) > 1 else ""
    return out


def load_consensus_meta(paths):
    # consensus fasta headers are stamped by RACON/MEDAKA as:
    #   >{sample}_{cluster_id} sample={sample} cluster_size={n_reads}
    rows = []
    for p in paths:
        with open(p) as fh:
            header = fh.readline().strip()
        if not header.startswith(">"):
            continue
        tokens = header[1:].split()
        seq_id = tokens[0]
        fields = dict(tok.split("=", 1) for tok in tokens[1:] if "=" in tok)
        rows.append({
            "path": p,
            "seq_id": seq_id,
            "sample": fields.get("sample"),
            "cluster_size": int(fields["cluster_size"]) if "cluster_size" in fields else None,
        })
    return pd.DataFrame(rows, columns=["path", "seq_id", "sample", "cluster_size"])


# flag_reason -> the consensus_by_confidence/ folder a cluster's fasta goes in
CONFIDENCE_DIR = {"": "confident", "low_identity": "low_confidence", "no_hit": "no_hit"}


def sort_consensus(best, out_dir):
    """Copy every consensus fasta into out_dir/{confident,low_confidence,no_hit}/
    as {sample}.{cluster_id}.consensus.fasta, by the same call the abundance
    table makes (flag_reason), so the two can't disagree."""
    for d in CONFIDENCE_DIR.values():
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    for row in best.dropna(subset=["path", "sample"]).itertuples():
        cluster_id = row.seq_id[len(row.sample) + 1:]  # seq_id is {sample}_{cluster_id}
        dest = os.path.join(out_dir, CONFIDENCE_DIR[row.flag_reason], f"{row.sample}.{cluster_id}.consensus.fasta")
        shutil.copyfile(row.path, dest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", nargs="+", required=True)
    ap.add_argument("--consensus", nargs="+", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-sorted", help="directory to gather consensus fastas into, by confidence")
    ap.add_argument("--min-pident", type=float, default=90,
                     help="hits below this identity are flagged low_identity, not dropped")
    ap.add_argument("--min-abundance", type=int, default=0,
                     help="merged clusters with fewer reads than this get low_abundance=true, not dropped")
    ap.add_argument("--min-rel-abundance", type=float, default=0,
                     help="merged clusters holding a smaller fraction of their sample's clustered reads "
                          "than this get low_abundance=true, not dropped")
    ap.add_argument("--silva-taxmap",
                     help="SILVA tax_slv_ssu file, to give each level of a SILVA lineage its true rank")
    args = ap.parse_args()

    taxmap = load_silva_taxmap(args.silva_taxmap) if args.silva_taxmap else None
    hits_df = load_hits(args.hits)
    consensus_meta = load_consensus_meta(args.consensus)

    # keep best hit per cluster (highest bitscore) as the working call.
    # Ties are common (several reference records with identical bitscores),
    # so the sort must be stable: BLAST lists a cluster's hits best-first, and
    # a stable sort keeps that order among tied bitscores, so the hit picked
    # here is the first line of the cluster's hits file. pandas' default sort
    # is not stable and picked an arbitrary tied hit, which gave clusters with
    # identical hits different species.
    best = (
        hits_df.sort_values("bitscore", ascending=False, kind="stable")
        .groupby("seq_id", as_index=False)
        .first()
    )
    best = best.merge(consensus_meta, on="seq_id", how="left")
    best["tied_taxa"] = best["seq_id"].map(tied_species(hits_df, taxmap)).fillna("")

    # lineage columns only when the reference is SILVA-style, so a RefSeq or
    # plain custom run's table keeps exactly its old columns
    lineage_cols = LINEAGE_COLS if best["stitle"].map(is_silva).any() else []
    if lineage_cols:
        lineages = best["stitle"].map(lambda s: silva_lineage(s, taxmap) if is_silva(s) else {})
        for col in lineage_cols:
            best[col] = lineages.map(lambda lin: lin.get(col, ""))

    # clusters are flagged, not dropped, either for no BLAST hit at all or for
    # a best hit too divergent to call with confidence (pident < min-pident)
    best["flag_reason"] = ""
    best.loc[best["subject_id"] == "NO_HIT", "flag_reason"] = "no_hit"
    low_pident = best["pident"].notna() & (best["pident"] < args.min_pident)
    best.loc[low_pident & (best["flag_reason"] == ""), "flag_reason"] = "low_identity"

    # Thin support is a separate axis from the BLAST call above, so it gets its
    # own true/false column rather than another flag_reason value: flag_reason
    # stays exactly '' | 'no_hit' | 'low_identity', which downstream readers
    # (e.g. the web frontend) take as the whole story and use as a taxon label.
    # Judged on the *merged* cluster size (reads from every cluster folded into
    # this one), so an organism split into several small read clusters isn't
    # penalised for the split. The relative form is a fraction of the sample's
    # clustered reads, which keeps the cutoff meaningful across barcodes of
    # very different depth.
    sample_total = best.groupby("sample")["cluster_size"].transform("sum")
    is_low_abundance = best["cluster_size"].notna() & (
        (best["cluster_size"] < args.min_abundance)
        | (best["cluster_size"] / sample_total < args.min_rel_abundance)
    )
    best["low_abundance"] = is_low_abundance.map({True: "true", False: "false"})

    if args.out_sorted:
        sort_consensus(best, args.out_sorted)

    # lineage columns sit where edna-ont-nf puts its --taxdump ones
    best = best[["seq_id", "sample", "cluster_size", "subject_id", "pident",
                 "length", "evalue", "bitscore", "stitle", "tied_taxa"] + lineage_cols
                + ["flag_reason", "low_abundance"]]
    best.to_csv(args.out_table, sep="\t", index=False)

    # "flagged" is the BLAST-confidence flag; low_abundance is counted
    # separately since it says how much support a call has, not whether the
    # call itself is doubtful
    n_clusters = best.shape[0]
    n_flagged = (best["flag_reason"] != "").sum()
    n_low_abundance = (best["low_abundance"] == "true").sum()
    per_sample = (
        best.groupby("sample")
        .agg(clusters=("seq_id", "count"),
             flagged=("flag_reason", lambda s: (s != "").sum()),
             low_abundance=("low_abundance", lambda s: (s == "true").sum()))
        .reset_index()
    )

    with open(args.out_html, "w") as fh:
        fh.write("<html><body><h2>Run summary</h2>")
        fh.write(f"<p>Clusters processed: {n_clusters}</p>")
        fh.write(f"<p>Clusters flagged for manual review (no hit or pident &lt; {args.min_pident}): {n_flagged}</p>")
        cutoffs = []
        if args.min_abundance:
            cutoffs.append(f"&lt; {args.min_abundance} reads")
        if args.min_rel_abundance:
            cutoffs.append(f"&lt; {args.min_rel_abundance:g} of the sample's clustered reads")
        if cutoffs:
            fh.write(f"<p>Clusters marked low_abundance ({' or '.join(cutoffs)}): {n_low_abundance}</p>")
        fh.write("<h3>Per-sample</h3>")
        fh.write(per_sample.to_html(index=False))
        fh.write("</body></html>")


if __name__ == "__main__":
    main()
