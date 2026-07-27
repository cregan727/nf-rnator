#!/usr/bin/env python3
"""
Convert a plate map's barcode column into a FASTA file for cutadapt-based
demultiplexing, using sample_name as the FASTA header (what cutadapt's
{name} template will use to name output files).

Standard-library only (no pandas) -- this runs inside the cutadapt
container alongside cutadapt itself, which won't have pandas installed.

Usage:
    build_barcode_fasta.py --platemap platemap_Library.csv --output barcodes.fasta
"""
import argparse
import csv
import re
import sys
from collections import Counter


def sanitize(name):
    """Keep filenames/FASTA headers safe: alphanumeric, dash, underscore, dot only."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--platemap", required=True, help="CSV: barcode,well,project,sample_name,genome")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.platemap, newline="") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        sys.exit(f"platemap {args.platemap} has no data rows")
    missing = {"barcode", "sample_name"} - set(rows[0].keys())
    if missing:
        sys.exit(f"platemap {args.platemap} is missing required column(s): {missing}")

    names = [sanitize(r["sample_name"]) for r in rows]
    counts = Counter(names)
    dupes = [n for n, c in counts.items() if c > 1]
    if dupes:
        sys.exit(f"platemap {args.platemap} has duplicate sample_name(s) after sanitizing "
                  f"for filename-safety: {dupes} -- cutadapt needs unique names per barcode.")

    with open(args.output, "w") as fh:
        for name, row in zip(names, rows):
            fh.write(f">{name}\n{row['barcode']}\n")

    print(f"Wrote {args.output}: {len(rows)} barcodes")


if __name__ == "__main__":
    main()
