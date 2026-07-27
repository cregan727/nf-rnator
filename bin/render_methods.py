#!/usr/bin/env python3
"""
Assemble a Methods section for one project: which genome(s)/annotation(s)
were used, the exact STARsolo command(s) that ran (from real captured
provenance records, not retyped), and a plain-language description of the
protocol (BRB-seq barcode/UMI structure, cell-filtering approach).

Only includes provenance for the (library, genome) combinations that
actually appear in this project's h5ad -- not every provenance file in the
run, in case other projects/plates used different genomes.

Usage:
    render_methods.py --h5ad ProjectX.h5ad --project ProjectX \
        --provenance LibA.mouse.provenance.txt LibA.human.provenance.txt ... \
        --output ProjectX_methods.txt
"""
import argparse
from pathlib import Path

import anndata as ad

BRBSEQ_DESCRIPTION = """\
Library preparation: Alithea MERCURIUS BRB-seq (bulk RNA barcoding and
sequencing), barcode set V5D. Read 1 carries a 14 bp cell/sample barcode
(positions 1-14) followed by a 14 bp UMI (positions 15-28); Read 2 carries
the cDNA sequence.

Alignment and quantification: reads were aligned and quantified with
STARsolo (--soloType CB_UMI_Simple), matching each read's 14 bp barcode to
the known set of sample barcodes with up to 1 mismatch
(--soloCBmatchWLtype 1MM), deduplicating UMIs with a directional
1-mismatch-aware algorithm (--soloUMIdedup 1MM_Directional). Unlike
droplet-based single-cell data, every barcode here corresponds to a known,
real bulk RNA sample rather than an unknown mixture of real cells and
empty droplets -- so no ambient-RNA/empty-droplet filtering was applied;
all barcodes present in the pool were retained (--soloCellFilter TopCells
96, --soloCellFilter None would achieve the same result on kits using
fewer than 96 wells).\
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True, help="Merged per-project .h5ad")
    ap.add_argument("--project", required=True)
    ap.add_argument("--provenance", required=True, nargs="+",
                     help="All *.provenance.txt files collected across the run -- "
                          "filtered down to just the ones relevant to this project")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    relevant_libraries = set(adata.obs["library"].unique()) if "library" in adata.obs else set()
    relevant_genomes = set(adata.obs["genome"].unique()) if "genome" in adata.obs else set()

    # provenance filenames are {library}.{genome}.provenance.txt
    relevant_files = []
    for f in args.provenance:
        stem = Path(f).name.removesuffix(".provenance.txt")
        parts = stem.rsplit(".", 1)
        if len(parts) != 2:
            continue
        library, genome = parts
        if library in relevant_libraries and genome in relevant_genomes:
            relevant_files.append(f)

    sections = [
        f"Methods -- {args.project}",
        "=" * (11 + len(args.project)),
        "",
        BRBSEQ_DESCRIPTION,
        "",
        f"This project's samples were processed on {len(relevant_files)} "
        f"plate/genome alignment run(s):",
        "",
    ]
    for f in sorted(relevant_files):
        sections.append(f"--- {Path(f).name} ---")
        sections.append(Path(f).read_text().rstrip())
        sections.append("")

    with open(args.output, "w") as fh:
        fh.write("\n".join(sections) + "\n")

    print(f"Wrote {args.output}: {len(relevant_files)} of {len(args.provenance)} "
          f"provenance record(s) were relevant to this project")


if __name__ == "__main__":
    main()
