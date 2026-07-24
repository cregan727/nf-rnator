#!/usr/bin/env python3
"""
Export a project-level .h5ad to a plain gene expression table: gene
SYMBOLS as rows, sample NAMES as columns, raw counts as values. This is
the standard format expected by several downstream tools (e.g. BioJupies'
upload spec), not a format specific to any one of them.

Also writes a companion metadata table (one row per sample: sample_name,
project, library, well, genome).

Usage:
    export_gene_table.py --h5ad ProjectX.h5ad --project ProjectX \
                         --counts-output ProjectX_counts.tsv \
                         --metadata-output ProjectX_metadata.tsv
"""
import argparse
import sys

import anndata as ad
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True, help="Merged per-project .h5ad from merge_project_h5ad.py")
    ap.add_argument("--project", required=True)
    ap.add_argument("--counts-output", required=True)
    ap.add_argument("--metadata-output", required=True)
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)

    if "gene_symbols" not in adata.var.columns:
        sys.exit(f"{args.h5ad} has no var['gene_symbols'] -- expected this to be populated "
                  f"automatically by scanpy's read_10x_mtx() from STARsolo's features.tsv. "
                  f"Check that the reference GTF used for alignment has gene_name attributes.")

    if "sample_name" not in adata.obs.columns:
        sys.exit(f"{args.h5ad} has no obs['sample_name'] -- expected this from the plate map.")

    # gene symbols as row labels -- fall back to the Ensembl ID (current
    # var_names) for any gene with a missing/blank symbol, then de-duplicate
    # (distinct genes can legitimately share a symbol, or a symbol can be
    # blank for several genes) since BioJupies expects unique row labels.
    symbols = adata.var["gene_symbols"].copy()
    blank = symbols.isna() | (symbols.astype(str).str.strip() == "")
    symbols[blank] = adata.var_names[blank]
    n_fallback = int(blank.sum())

    counts = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X),
        index=adata.obs["sample_name"].values,
        columns=symbols.values,
    ).T  # genes as rows, samples as columns, per BioJupies' spec
    counts.index.name = "gene_symbol"
    # raw counts are whole numbers even though anndata stores them as float;
    # round rather than truncate in case of any upstream floating-point noise
    counts = counts.round().astype(int)

    # de-duplicate row labels (keep all rows, just make labels unique --
    # dropping/summing duplicate-symbol genes would silently lose data)
    if counts.index.duplicated().any():
        counts.index = pd.io.common.dedup_names(counts.index, is_potential_multiindex=False) \
            if hasattr(pd.io.common, "dedup_names") else _dedup(counts.index)

    counts.to_csv(args.counts_output, sep="\t")

    meta_cols = [c for c in ["sample_name", "project", "library", "well", "genome", "n_reads", "n_counts"]
                 if c in adata.obs.columns]
    metadata = adata.obs[meta_cols].drop_duplicates(subset="sample_name").set_index("sample_name")
    metadata.to_csv(args.metadata_output, sep="\t")

    print(f"Wrote {args.counts_output}: {counts.shape[0]} genes x {counts.shape[1]} samples "
          f"({n_fallback} gene(s) fell back to Ensembl ID, no symbol available)")
    print(f"Wrote {args.metadata_output}: {len(metadata)} samples")


def _dedup(index):
    """Fallback de-duplication for older pandas without io.common.dedup_names."""
    seen = {}
    out = []
    for name in index:
        if name in seen:
            seen[name] += 1
            out.append(f"{name}.{seen[name]}")
        else:
            seen[name] = 0
            out.append(name)
    return pd.Index(out)


if __name__ == "__main__":
    main()
