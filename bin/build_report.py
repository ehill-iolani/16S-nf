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
import pandas as pd

COLS = ["seq_id", "subject_id", "pident", "length", "evalue", "bitscore", "stitle"]


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


def species_label(stitle):
    # same rule as the frontend's shortTaxonLabel: BLAST's stitle starts with
    # the binomial, so the first two tokens name the species
    toks = stitle.split() if isinstance(stitle, str) else []
    return " ".join(toks[:2])


def tied_species(hits_df):
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
        names = list(dict.fromkeys(n for n in map(species_label, top["stitle"]) if n))
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
            "seq_id": seq_id,
            "sample": fields.get("sample"),
            "cluster_size": int(fields["cluster_size"]) if "cluster_size" in fields else None,
        })
    return pd.DataFrame(rows, columns=["seq_id", "sample", "cluster_size"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", nargs="+", required=True)
    ap.add_argument("--consensus", nargs="+", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--min-pident", type=float, default=90,
                     help="hits below this identity are flagged low_identity, not dropped")
    ap.add_argument("--min-abundance", type=int, default=0,
                     help="merged clusters with fewer reads than this get low_abundance=true, not dropped")
    ap.add_argument("--min-rel-abundance", type=float, default=0,
                     help="merged clusters holding a smaller fraction of their sample's clustered reads "
                          "than this get low_abundance=true, not dropped")
    args = ap.parse_args()

    hits_df = load_hits(args.hits)
    consensus_meta = load_consensus_meta(args.consensus)

    # keep best hit per cluster (highest bitscore) as the working call.
    # Ties are common (several reference records with identical bitscores),
    # so the sort must be stable: BLAST lists a cluster's hits best-first, and
    # a stable sort keeps that order among tied bitscores, so the hit picked
    # here is the first line of the cluster's hits file -- the same one
    # SORT_CONSENSUS reads to choose the confidence folder. pandas' default
    # sort is not stable and picked an arbitrary tied hit, which could disagree
    # with SORT_CONSENSUS (a cluster flagged low_identity but filed under
    # confident/) and gave clusters with identical hits different species.
    best = (
        hits_df.sort_values("bitscore", ascending=False, kind="stable")
        .groupby("seq_id", as_index=False)
        .first()
    )
    best = best.merge(consensus_meta, on="seq_id", how="left")
    best["tied_taxa"] = best["seq_id"].map(tied_species(hits_df)).fillna("")

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

    best = best[["seq_id", "sample", "cluster_size", "subject_id", "pident",
                 "length", "evalue", "bitscore", "stitle", "tied_taxa", "flag_reason", "low_abundance"]]
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
