#!/usr/bin/env python3
"""
Merge one plate's per-genome h5ads into a single whole-plate h5ad.

Each fragment (one per genome present on the plate) already had every well
NOT assigned to that genome dropped by mtx_to_h5ad.py, so the fragments
cover disjoint sets of wells (and disjoint blocks of genes, since different
species) -- concatenating them just unions those cleanly, with zero-padding
for wells/genes that don't apply to a given fragment.

Usage:
    merge_plate_genomes.py --inputs LibraryA.mouse.h5ad LibraryA.human.h5ad \
                            --library LibraryA --output LibraryA.h5ad
"""
import argparse
import anndata as ad


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", required=True, nargs="+", help="Per-genome h5ad fragments for this plate")
    ap.add_argument("--library", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if len(args.inputs) == 1:
        adata = ad.read_h5ad(args.inputs[0])
    else:
        fragments = [ad.read_h5ad(f) for f in args.inputs]
        adata = ad.concat(fragments, join="outer", merge="first", index_unique=None)

    adata.write_h5ad(args.output)
    genomes = sorted(adata.obs["genome"].unique()) if "genome" in adata.obs else ["?"]
    print(f"Wrote {args.output}: {adata.n_obs} wells x {adata.n_vars} genes "
          f"from {len(args.inputs)} genome(s) {genomes} (library={args.library})")


if __name__ == "__main__":
    main()
