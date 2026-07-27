#!/usr/bin/env python3
"""
Move cutadapt's flat per-sample fastq output into project subdirectories,
based on the plate map's sample_name -> project mapping. Filenames get
prefixed with the library name, so the same sample_name reused across
different plates/libraries can't collide.

Standard-library only (no pandas) -- this runs inside the cutadapt
container alongside cutadapt itself, which won't have pandas installed.

Usage:
    route_by_project.py --platemap platemap_Library.csv --library Library --indir . --outdir projects
"""
import argparse
import csv
import re
import shutil
import sys
from pathlib import Path


def sanitize(name):
    """Must match build_barcode_fasta.py's sanitizing exactly, or filenames won't line up."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--platemap", required=True)
    ap.add_argument("--library", required=True)
    ap.add_argument("--indir", required=True, help="Directory with cutadapt's flat {sample_name}_R1/R2.fastq.gz output")
    ap.add_argument("--outdir", required=True, help="Root directory to organize into {project}/fastqs/")
    args = ap.parse_args()

    with open(args.platemap, newline="") as fh:
        rows = list(csv.DictReader(fh))

    indir = Path(args.indir)
    outdir = Path(args.outdir)

    moved = 0
    for row in rows:
        name = sanitize(row["sample_name"])
        project = sanitize(row["project"])
        for read in ("R1", "R2"):
            src = indir / f"{name}_{read}.fastq.gz"
            if not src.exists():
                print(f"WARNING: expected output {src} not found (0 reads for this sample?)", file=sys.stderr)
                continue
            dest_dir = outdir / project / "fastqs"
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest = dest_dir / f"{args.library}_{name}_{read}.fastq.gz"
            shutil.move(str(src), str(dest))
            moved += 1

    # keep the "unknown" (unmatched-barcode) bucket too, for QC visibility --
    # not attributed to any project since it isn't attributable to any sample
    for read in ("R1", "R2"):
        src = indir / f"unknown_{read}.fastq.gz"
        if src.exists():
            dest_dir = outdir / "_unmatched"
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dest_dir / f"{args.library}_unknown_{read}.fastq.gz"))

    print(f"Moved {moved} fastq files into project subdirectories under {outdir}")


if __name__ == "__main__":
    main()
