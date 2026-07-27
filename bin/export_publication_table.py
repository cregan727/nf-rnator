#!/usr/bin/env python3
"""
Build a publication-ready gene expression table: versioned Ensembl ID as
the row key, one column per sample (raw counts), plus the full GENCODE
gene-annotation columns (seqnames, start, end, width, strand, source,
gene_type, gene_name, hgnc_id, havana_gene) -- matching the standard
format requested by journals/reviewers (a GENCODE-annotated SummarizedExperiment
rowData table, essentially).

Genes in the h5ad with no matching row in the annotation table (e.g. a
project spans multiple genomes, but this table only has one genome's
annotation) get blank annotation columns rather than being silently
dropped -- reported as a count so it's not silently lossy.

Usage:
    export_publication_table.py --h5ad ProjectX.h5ad \
        --gene-annotation mouse_annotation.tsv human_annotation.tsv \
        --project ProjectX --output ProjectX_publication_table.tsv
"""
import argparse
import csv
import sys

import anndata as ad
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True, help="Merged per-project .h5ad")
    ap.add_argument("--gene-annotation", required=True, nargs="+",
                     help="gene_annotation.tsv from gene_annotation_table.py -- one per genome "
                          "present in this project (Ensembl IDs are globally unique across "
                          "species, so these just get combined into one lookup)")
    ap.add_argument("--project", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    if "sample_name" not in adata.obs.columns:
        sys.exit(f"{args.h5ad} has no obs['sample_name'] -- expected this from the plate map.")

    # combine every genome's annotation table into one lookup (Ensembl IDs
    # are globally unique across species, so no collision risk)
    annotation = {}
    annot_fields = None
    for path in args.gene_annotation:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if annot_fields is None:
                annot_fields = [f for f in reader.fieldnames if f != "EnsemblID"]
            for row in reader:
                annotation[row["EnsemblID"]] = row

    sample_names = adata.obs["sample_name"].tolist()
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X = np.rint(X).astype(int)  # raw counts are whole numbers

    fieldnames = ["EnsemblID"] + sample_names + annot_fields
    n_missing = 0

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene_idx, gene_id in enumerate(adata.var_names):
            row = {"EnsemblID": gene_id}
            for sample_idx, sample in enumerate(sample_names):
                row[sample] = int(X[sample_idx, gene_idx])
            annot_row = annotation.get(gene_id)
            if annot_row is None:
                n_missing += 1
                for f in annot_fields:
                    row[f] = ""
            else:
                for f in annot_fields:
                    row[f] = annot_row[f]
            writer.writerow(row)

    print(f"Wrote {args.output}: {adata.n_vars} genes x {len(sample_names)} samples "
          f"({n_missing} gene(s) had no matching annotation row -- left blank, not dropped)")


if __name__ == "__main__":
    main()
