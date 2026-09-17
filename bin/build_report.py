#!/usr/bin/env python3
"""
Aggregate per-cluster BLAST hits + consensus fastas into a single
sample x taxon abundance table, plus a minimal QC/summary HTML.

This is a functional starting point, not the final report -- extend with
per-cluster read counts (from isONclust final_clusters.tsv), top-hit
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
    args = ap.parse_args()

    hits_df = load_hits(args.hits)
    consensus_meta = load_consensus_meta(args.consensus)

    # keep best hit per cluster (highest bitscore) as the working call
    best = (
        hits_df.sort_values("bitscore", ascending=False)
        .groupby("seq_id", as_index=False)
        .first()
    )
    best = best.merge(consensus_meta, on="seq_id", how="left")

    # clusters are flagged, not dropped, either for no BLAST hit at all or for
    # a best hit too divergent to call with confidence (pident < min-pident)
    best["flag_reason"] = ""
    best.loc[best["subject_id"] == "NO_HIT", "flag_reason"] = "no_hit"
    low_pident = best["pident"].notna() & (best["pident"] < args.min_pident)
    best.loc[low_pident & (best["flag_reason"] == ""), "flag_reason"] = "low_identity"

    best = best[["seq_id", "sample", "cluster_size", "subject_id", "pident",
                 "length", "evalue", "bitscore", "stitle", "flag_reason"]]
    best.to_csv(args.out_table, sep="\t", index=False)

    n_clusters = best.shape[0]
    n_flagged = (best["flag_reason"] != "").sum()
    per_sample = (
        best.groupby("sample")
        .agg(clusters=("seq_id", "count"), flagged=("flag_reason", lambda s: (s != "").sum()))
        .reset_index()
    )

    with open(args.out_html, "w") as fh:
        fh.write("<html><body><h2>Run summary</h2>")
        fh.write(f"<p>Clusters processed: {n_clusters}</p>")
        fh.write(f"<p>Clusters flagged for manual review (no hit or pident &lt; {args.min_pident}): {n_flagged}</p>")
        fh.write("<h3>Per-sample</h3>")
        fh.write(per_sample.to_html(index=False))
        fh.write("</body></html>")


if __name__ == "__main__":
    main()
