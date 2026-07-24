#!/usr/bin/env python3
"""
Convert a STARsolo raw count matrix to .h5ad, annotating each barcode with
its project and sample name from the library's plate map.

If a plate has multiple genomes present, this gets run once per genome
(against that genome's own STARsolo alignment of the WHOLE plate) with
--genome set -- only wells assigned to THAT genome in the plate map are
kept in the .h5ad; every other well is dropped there, since gene-level data
for them is meaningless (wrong-species alignment). The per-genome .h5ad
fragments get recombined into one correct whole-plate h5ad downstream by
merge_plate_genomes.py.

Separately, a lightweight --well-totals-output CSV is written BEFORE that
filtering happens, capturing total UMI count for EVERY well (not just
assigned ones) from THIS genome's alignment. That's what lets the plate
report show, per well, how many UMIs came from its assigned genome vs. any
other genome present on the plate -- a well with unexpectedly high
off-target-genome counts is a real mixup/contamination signal that would
otherwise be silently discarded.

--cell-reads-stats takes STARsolo's CellReads.stats (from
--soloCellReadStats Standard), which carries real RAW read counts per
barcode (the "cbMatch" column: reads whose CB matched this barcode) --
distinct from n_counts (deduplicated UMI count), which is all that was
available before. Optional: if not given, n_reads just won't be populated.

Usage:
    mtx_to_h5ad.py --solo-dir Library.mouse.Solo.out/Gene/raw \
                    --cell-reads-stats Library.mouse.Solo.out/Gene/CellReads.stats \
                    --platemap platemap_Library.csv \
                    --library Library --genome mouse \
                    --output Library.mouse.h5ad \
                    --well-totals-output Library.mouse.well_totals.csv
"""
import argparse
import sys
import numpy as np
import scanpy as sc
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solo-dir", required=True, help="STARsolo Solo.out/Gene/raw directory")
    ap.add_argument("--cell-reads-stats", default=None,
                     help="STARsolo Solo.out/Gene/CellReads.stats (from --soloCellReadStats Standard) -- "
                          "optional, adds real raw read counts (n_reads) per well")
    ap.add_argument("--platemap", required=True, help="CSV: barcode,well,project,sample_name,genome")
    ap.add_argument("--library", required=True, help="Library/plate name")
    ap.add_argument("--genome", required=True, help="Genome key this alignment was run against")
    ap.add_argument("--output", required=True, help="Output .h5ad path")
    ap.add_argument("--well-totals-output", required=True,
                     help="Output CSV: total UMIs per well from THIS genome run, before filtering")
    args = ap.parse_args()

    adata = sc.read_10x_mtx(args.solo_dir, var_names="gene_ids")
    adata.var_names_make_unique()

    # real raw read counts per barcode (pre-deduplication) -- distinct from
    # n_counts, which is deduplicated UMI count
    read_counts = None
    if args.cell_reads_stats:
        try:
            stats = pd.read_csv(args.cell_reads_stats, sep=r"\s+")
            read_counts = stats.set_index("CB")["cbMatch"]
        except Exception as e:
            print(f"WARNING: couldn't parse {args.cell_reads_stats} ({e}) -- "
                  f"n_reads will be unavailable", file=sys.stderr)

    platemap_full = pd.read_csv(args.platemap, dtype=str)
    required_cols = {"barcode", "well", "project", "sample_name", "genome"}
    missing = required_cols - set(platemap_full.columns)
    if missing:
        sys.exit(f"platemap {args.platemap} is missing required column(s): {missing}")

    # --- well_totals: every well present in this genome's raw matrix, ---
    # --- regardless of which genome it's actually assigned to ---
    total_umis = np.asarray(adata.X.sum(axis=1)).ravel()
    totals = pd.DataFrame({"barcode": adata.obs_names, "umis_this_genome": total_umis})
    totals = totals.merge(
        platemap_full.rename(columns={"genome": "assigned_genome"})[
            ["barcode", "well", "project", "sample_name", "assigned_genome"]
        ],
        on="barcode", how="left",
    )
    totals["this_genome"] = args.genome
    totals.to_csv(args.well_totals_output, index=False)

    # --- filtered h5ad: only wells actually assigned to this genome ---
    platemap = platemap_full[platemap_full["genome"] == args.genome].set_index("barcode")
    if len(platemap) == 0:
        sys.exit(f"No wells in {args.platemap} are assigned genome='{args.genome}' -- "
                  f"nothing for this alignment to keep.")

    adata.obs["library"] = args.library
    adata.obs["barcode"] = adata.obs_names
    adata.obs = adata.obs.join(platemap[["well", "project", "sample_name"]], on="barcode")
    if read_counts is not None:
        adata.obs["n_reads"] = adata.obs["barcode"].map(read_counts)

    # drop every well NOT assigned to this genome (expected/normal on a
    # mixed-genome plate -- those wells' real gene-level data comes from a
    # different genome's alignment of this same plate; their cross-genome
    # UMI totals were already captured above for QC purposes)
    n_before = adata.n_obs
    adata = adata[~adata.obs["project"].isna()].copy()
    print(f"{args.library}/{args.genome}: kept {adata.n_obs} of {n_before} wells "
          f"(rest belong to other genomes on this plate)", file=sys.stderr)

    adata.obs["genome"] = args.genome
    # obs_names unique across libraries AND genomes: same barcode sequence
    # can recur in a different plate, or even the same plate under a
    # different genome key before merging drops the wrong-genome copies.
    adata.obs_names = args.library + "_" + args.genome + "_" + adata.obs_names

    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}: {adata.n_obs} barcodes x {adata.n_vars} genes "
          f"({adata.obs['project'].nunique(dropna=True)} project(s))")


if __name__ == "__main__":
    main()
