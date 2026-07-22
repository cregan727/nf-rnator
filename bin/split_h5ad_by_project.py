#!/usr/bin/env python3
"""
Split a library-level (per-plate) .h5ad into one .h5ad per project, based on
the obs['project'] annotation set by mtx_to_h5ad.py.

Usage:
    split_h5ad_by_project.py --input Library.h5ad --library Library --outdir split/
"""
import argparse
import sys
from pathlib import Path
import anndata as ad


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Library-level .h5ad from mtx_to_h5ad.py")
    ap.add_argument("--library", required=True, help="Library/plate name (used in output filenames)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    adata = ad.read_h5ad(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    projects = adata.obs["project"].dropna().unique()
    if len(projects) == 0:
        sys.exit(f"No barcodes in {args.input} have a project assigned -- check the plate map.")

    for project in projects:
        subset = adata[adata.obs["project"] == project].copy()
        # sanitize project name for use in a filename
        safe_project = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(project))
        out_path = outdir / f"{args.library}__{safe_project}.h5ad"
        subset.write_h5ad(out_path)
        print(f"Wrote {out_path}: {subset.n_obs} barcodes (project={project})")


if __name__ == "__main__":
    main()
