#!/usr/bin/env python3
"""
Merge one project's h5ad fragments (potentially one per library/plate, if
that project's samples were spread across multiple sequencing runs) into a
single final .h5ad for that project.

Usage:
    merge_project_h5ad.py --inputs frag1.h5ad frag2.h5ad --project ProjectX --output ProjectX.h5ad
"""
import argparse
import anndata as ad


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", required=True, nargs="+", help="Per-library h5ad fragments for this project")
    ap.add_argument("--project", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if len(args.inputs) == 1:
        adata = ad.read_h5ad(args.inputs[0])
    else:
        fragments = [ad.read_h5ad(f) for f in args.inputs]
        # outer join: if this project's samples span libraries aligned to
        # DIFFERENT genomes, gene sets won't match -- outer join keeps all
        # genes, padding with zeros where a gene doesn't exist in a given
        # library's reference. Mixed-genome projects are a real edge case
        # here; the resulting matrix will have disjoint blocks of
        # mouse/human genes rather than anything biologically merged.
        adata = ad.concat(fragments, join="outer", label="library_batch", index_unique=None)

    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}: {adata.n_obs} barcodes x {adata.n_vars} genes "
          f"from {len(args.inputs)} librar{'y' if len(args.inputs)==1 else 'ies'} (project={args.project})")


if __name__ == "__main__":
    main()
