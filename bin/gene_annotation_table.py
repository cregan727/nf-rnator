#!/usr/bin/env python3
"""
Parse a GENCODE-format reference GTF into a gene-level annotation table:
one row per gene, versioned Ensembl ID as the key, plus seqnames, start,
end, width, strand, source, gene_type, gene_name, hgnc_id, havana_gene --
the standard GENCODE gene-annotation columns (matches what
rtracklayer::import() + subset(type=="gene") gives you in R/Bioconductor).

This information is NOT recoverable from a STAR genome index or a BAM file
-- STAR's own index only keeps gene_id/gene_name/biotype internally, none
of the rest. It has to come from the original reference GTF (the same file
the STAR index was built from).

hgnc_id and havana_gene are only present on HAVANA-source entries in real
GENCODE GTFs (ENSEMBL-source entries lack them) -- those columns are left
blank for genes that don't have them, not dropped.

Usage:
    gene_annotation_table.py --gtf annotation.gtf --output gene_annotation.tsv
"""
import argparse
import csv
import gzip
import re
import sys

ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_attributes(field):
    return dict(ATTR_RE.findall(field))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="Reference GTF (GENCODE format), .gtf or .gtf.gz")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    fieldnames = ["EnsemblID", "seqnames", "start", "end", "width", "strand",
                  "source", "gene_type", "gene_name", "hgnc_id", "havana_gene"]
    n_genes = 0
    n_missing_hgnc = 0

    with open_maybe_gz(args.gtf) as fh, open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            seqname, source, _, start, end, _, strand, _, attr_field = fields
            attrs = parse_attributes(attr_field)
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            start_i, end_i = int(start), int(end)
            row = {
                "EnsemblID": gene_id,
                "seqnames": seqname,
                "start": start_i,
                "end": end_i,
                "width": end_i - start_i + 1,
                "strand": strand,
                "source": source,
                "gene_type": attrs.get("gene_type", attrs.get("gene_biotype", "")),
                "gene_name": attrs.get("gene_name", ""),
                "hgnc_id": attrs.get("hgnc_id", ""),
                "havana_gene": attrs.get("havana_gene", ""),
            }
            if not row["hgnc_id"]:
                n_missing_hgnc += 1
            writer.writerow(row)
            n_genes += 1

    if n_genes == 0:
        sys.exit(f"No 'gene' feature lines found in {args.gtf} -- is this a valid GTF? "
                  f"(GFF3 uses different field 3 conventions and isn't handled here.)")

    print(f"Wrote {args.output}: {n_genes} genes "
          f"({n_missing_hgnc} without an hgnc_id -- expected for non-HAVANA-source entries)")


if __name__ == "__main__":
    main()
