#!/usr/bin/env python3
"""
Bray-Curtis dissimilarity + classical PCoA across samples, computed from the
same abundance_table.tsv BUILD_REPORT produces. Implemented with pandas +
numpy only (no scipy/scikit-bio) so it runs in the same biocontainer
BUILD_REPORT already uses -- numpy comes along as pandas' own runtime
dependency, so nothing new has to be verified to work on Batch.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def short_taxon_label(stitle, subject_id):
    """Mirrors resultsFormat.js's shortTaxonLabel on the frontend -- both
    must agree on what counts as "one species" so a taxon in the abundance
    chart/rarefaction curve corresponds to the same grouping used here for
    beta diversity. Keep these two implementations in sync if either changes.
    """
    words = str(stitle).split() if pd.notna(stitle) else []
    if len(words) >= 2:
        return f"{words[0]} {words[1]}"
    if len(words) == 1:
        return words[0]
    return subject_id if pd.notna(subject_id) else "Unclassified"


def resolved_taxon(row):
    """The taxon a confidently-hit cluster counts as. Hits tied for the best
    bitscore that name several species (build_report.py's tied_taxa, "|"-joined)
    are called at genus level when they share a genus ("Genus sp.") and
    dropped (None) when they don't; the frontend applies the same rule."""
    tied = row.get("tied_taxa")
    names = tied.split("|") if isinstance(tied, str) and tied else []
    if len(names) > 1:
        genera = {name.split()[0] for name in names}
        return f"{genera.pop()} sp." if len(genera) == 1 else None
    return short_taxon_label(row.get("stitle"), row.get("subject_id"))


def bray_curtis(matrix):
    """matrix: samples x taxa abundance array. Returns an n x n distance matrix.
    Sample counts here are classroom-scale (tens, not thousands), so the
    plain O(n^2) loop is simpler than vectorizing and plenty fast."""
    n = matrix.shape[0]
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            denom = (matrix[i] + matrix[j]).sum()
            d = np.abs(matrix[i] - matrix[j]).sum() / denom if denom > 0 else 0.0
            dist[i, j] = dist[j, i] = d
    return dist


def classical_pcoa(dist, n_axes):
    """Gower's classical PCoA: double-center the squared distance matrix and
    eigendecompose. Negative eigenvalues (common with non-Euclidean distances
    like Bray-Curtis) are dropped rather than producing imaginary
    coordinates -- standard practice, not a bug."""
    n = dist.shape[0]
    d2 = dist**2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j

    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = np.clip(eigenvalues, 0, None)
    total = positive.sum()

    # A centered distance matrix for n samples has rank at most n-1, so
    # there's nothing meaningful past that many axes.
    usable_axes = min(n_axes, n - 1) if n > 1 else 0
    coords = np.zeros((n, usable_axes))
    pct_explained = np.zeros(usable_axes)
    for k in range(usable_axes):
        coords[:, k] = eigenvectors[:, k] * np.sqrt(positive[k])
        pct_explained[k] = (positive[k] / total * 100) if total > 0 else 0.0
    return coords, pct_explained


def empty_output(n_axes):
    cols = ["sample"] + [f"PC{i + 1}" for i in range(n_axes)] + [f"PC{i + 1}_pct" for i in range(n_axes)]
    return pd.DataFrame(columns=cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", required=True, help="abundance_table.tsv from BUILD_REPORT")
    ap.add_argument("--metadata", help="optional CSV: sample,<metadata columns>")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-axes", type=int, default=3)
    args = ap.parse_args()

    df = pd.read_csv(args.abundance, sep="\t")
    samples = sorted(df["sample"].dropna().unique())

    if len(samples) < 2:
        # PCoA needs at least two samples to place on axes -- write a
        # header-only file rather than fail the whole run over it.
        empty_output(args.n_axes).to_csv(args.out, sep="\t", index=False)
        return

    # Only confidently-identified taxa count toward community composition --
    # the same exclusion the frontend's rarefaction chart applies, so the two
    # analyses of the same run stay consistent with each other.
    resolved = df[df["flag_reason"].fillna("") == ""].copy()
    resolved["taxon"] = resolved.apply(resolved_taxon, axis=1)
    # clusters tied between species of different genera are "ambiguous" in the
    # frontend, so they don't count here either
    resolved = resolved[resolved["taxon"].notna()]

    matrix_df = (
        resolved.groupby(["sample", "taxon"])["cluster_size"]
        .sum()
        .unstack(fill_value=0)
        .reindex(samples, fill_value=0)
    )

    dist = bray_curtis(matrix_df.to_numpy(dtype=float))
    coords, pct_explained = classical_pcoa(dist, args.n_axes)

    out = pd.DataFrame({"sample": samples})
    for i in range(args.n_axes):
        if i < coords.shape[1]:
            out[f"PC{i + 1}"] = coords[:, i]
            out[f"PC{i + 1}_pct"] = pct_explained[i]
        else:
            # Fewer usable axes than requested (very small sample counts) --
            # pad with zeros so the output always has the same fixed set of
            # columns for the frontend to parse.
            out[f"PC{i + 1}"] = 0.0
            out[f"PC{i + 1}_pct"] = 0.0

    if args.metadata:
        meta = pd.read_csv(args.metadata, dtype=str)
        if "sample" not in meta.columns:
            sys.exit("--metadata CSV must have a 'sample' column")
        out = out.merge(meta, on="sample", how="left")

    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
