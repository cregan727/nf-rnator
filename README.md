# brbseq-pipeline

A minimal Nextflow pipeline for Alithea MERCURIUS BRB-seq data (kit PN 10813,
V5A barcode set), built to do only what BRB-seq actually needs rather than
carry the generality of a pipeline like nf-core/scrnaseq (many aligners, a
10x-oriented resource/schema system, ambient-RNA removal tuned for droplet
data). One aligner (STARsolo), one code path, full transparency into every
command that actually runs.

FastQC → multi-genome STARsolo (a plate can mix species across wells) →
per-well h5ad conversion → split/merge by project → per-plate and
per-project HTML QC reports with a real 8x12 spatial plate map, including
detection of cross-genome sample mixups.

## Requirements

- [Nextflow](https://www.nextflow.io/) (developed against 26.04.x)
- Singularity or Apptainer
- A prebuilt STAR genome index per genome you'll align against (this
  pipeline does not build indices for you -- see "Genome indices" below)

## Setup

```bash
git clone <this-repo>
cd brbseq-pipeline
chmod +x bin/*.py    # REQUIRED -- see note below
```

**Why `chmod +x` is required:** Nextflow calls everything in `bin/` by bare
filename (e.g. `mtx_to_h5ad.py --input ...`), which only works if those
files are executable and on `PATH` (Nextflow adds `bin/` to `PATH`
automatically, but doesn't set the executable bit for you). If this repo
reached you through a route that doesn't preserve file permissions --
downloaded from a chat interface, some non-git zip/tar packaging, etc. --
this bit will be lost even though it's set correctly in the git history, so
run this once after cloning/downloading regardless.

## Quick start

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genomes genomes.csv \
    --whitelist V5A_barcodes.txt \
    --outdir results
```

See `examples/` for the exact format of each input file, described below.

## Input files

### 1. `--input`: samplesheet.csv

One row per pooled library/plate. Genome is **not** here -- see plate map
below, since a single plate can mix species across wells.

```csv
library,fastq_1,fastq_2,plate_map
LibraryA,/path/to/LibraryA_R1.fastq.gz,/path/to/LibraryA_R2.fastq.gz,platemap_LibraryA.csv
LibraryB,/path/to/LibraryB_R1.fastq.gz,/path/to/LibraryB_R2.fastq.gz,platemap_LibraryB.csv
```

- `fastq_1` = R1 (14bp cell barcode + 14bp UMI), `fastq_2` = R2 (cDNA) --
  standard 10x-style read layout, which BRB-seq shares.
- If a library was sequenced across multiple lanes, merge the lane fastqs
  into one R1/R2 pair per library before running this (e.g. `cat
  *_L00*_R1_001.fastq.gz > Library_R1.fastq.gz`), or extend `main.nf` to
  accept comma-joined lane lists -- not currently built in.

### 2. `--genomes`: genomes.csv

Maps a short genome key to a prebuilt STAR index directory.

```csv
genome,star_index
mouse,/path/to/mouse/star_index
human,/path/to/human/star_index
```

### 3. Plate map (referenced per-library from the samplesheet)

One row per **used well**, mapping barcode sequence to physical well
position, project, sample name, and genome.

```csv
barcode,well,project,sample_name,genome
TACGTTATTCCGAA,A01,ProjectX,SampleX_01,mouse
AACAGGATAACTCC,B01,ProjectX,SampleX_02,mouse
ACTCAGGCACCTCC,C01,ProjectY,SampleY_01,human
```

- `well`: physical plate position, `A01`-`H12`. Required for the spatial
  plate-map visualizations in the report; if omitted, those sections are
  skipped gracefully.
- `genome`: can vary row-by-row -- a single plate mixing e.g. mouse and
  human wells is supported (see "Mixed-genome plates" below). If every well
  on a plate is the same species, just repeat that value on every row.
- A plate belonging entirely to one project is just every row sharing the
  same `project` value -- this format doesn't force per-well granularity,
  it just allows for it.
- You don't need every one of the 96 possible wells listed, only the ones
  actually used in that pool.

### 4. `--whitelist`: barcode whitelist

The full V5A barcode set for this kit -- one 14bp sequence per line, no
header. Get this from Alithea for your specific kit/lot (they distribute it
via a gated form/support email, not a public download). See
`examples/` for the expected plain format (not included here, since the
actual sequences are kit-specific).

## Genome indices

This pipeline expects a prebuilt STAR index per genome, referenced by path
in `genomes.csv`. Build one with:

```bash
STAR --runMode genomeGenerate --runThreadN 16 \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf
```

If you already have a 10x Cell Ranger reference package, don't reuse its
bundled `star/` index directly -- Cell Ranger pins its own STAR build
internally, and a version mismatch against whatever STAR you're running
here can make STAR reject the index outright. Build fresh from that
package's `fasta/genome.fa` + `genes/genes.gtf` (or equivalent) instead.

## Mixed-genome plates

STAR aligns a whole fastq file against one genome index per run -- it has
no concept of "use this index for these barcodes, that index for those."
So for any plate, this pipeline determines every distinct genome present
(from that plate's map) and runs STARsolo once per genome, aligning the
**whole plate's reads** against each one. `mtx_to_h5ad.py` then keeps only
the wells actually assigned to that run's genome and drops the rest;
`merge_plate_genomes.py` recombines the (disjoint) per-genome fragments
into one correct whole-plate h5ad.

**Cost tradeoff:** a plate with N genomes present costs roughly Nx the
compute of a single-genome plate -- every genome's STARsolo run processes
100% of that plate's reads, not a proportional subset. This trades compute
for avoiding a separate, fragile pre-alignment demultiplexing-by-raw-barcode
step. Fine for occasional mixed plates; worth knowing if mixed plates
become the norm.

**Mixup detection:** because every genome's run sees the whole plate before
filtering, the pipeline also captures each well's UMI count from every
genome, not just its assigned one. The plate report flags wells where more
than 10% of total UMIs came from an unassigned genome -- a real signal of
sample swap or cross-contamination -- while ignoring wells with too little
total signal (<500 UMIs) for that ratio to be statistically meaningful.

## Known limitations

- If a single project's samples are split across plates aligned to
  *different* genomes, the merged per-project h5ad will contain disjoint
  blocks of genes from each genome (outer-joined, zero-padded), not
  anything biologically reconciled across species.
- No downstream analysis (normalization, clustering, UMAP, differential
  expression) -- this pipeline stops at QC'd, correctly-attributed count
  matrices. BRB-seq samples are bulk RNA samples multiplexed via barcode,
  not individual cells, so single-cell-style clustering/UMAP isn't a
  meaningful next step anyway; you'd want a standard bulk differential
  expression workflow (DESeq2/edgeR-style) on the output matrices instead.
- Container image pins in `nextflow.config` (particularly the scanpy one)
  are a reasonable first guess, not independently verified against a live
  registry -- check they still resolve before relying on this at scale.
- Not run against real data end-to-end by the author of this scaffold --
  the Python components (`bin/*.py`) were tested against synthetic
  fixtures during development (including a deliberately-injected mixup
  scenario, confirmed correctly detected), but the full `main.nf` DAG
  itself has not been executed. Treat the first real run as a genuine test,
  not a guaranteed-correct deliverable -- run on one small plate before
  trusting it broadly.

## Output structure

```
results/
├── plates/
│   └── <library>/
│       ├── <library>.h5ad              # whole-plate, all genomes combined
│       ├── <library>_report.html       # per-plate QC report
│       ├── <genome>/                   # per-genome STARsolo output, logs, fastqc
│       └── fastqc/
└── projects/
    └── <project>/
        ├── <project>.h5ad              # merged across every plate contributing to this project
        └── <project>_report.html       # per-project QC report
```

## Running on SLURM

See `submit_slurm.sh` for a template. Key points baked into it from hands-on
debugging against a real cluster:

- `process.executor = 'local'` is set explicitly in `nextflow.config` so
  everything runs inside one job's allocation -- no Nextflow-submits-its-
  own-sub-jobs behavior.
- Path resolution uses `$SLURM_SUBMIT_DIR`, not `${BASH_SOURCE[0]}` --
  some SLURM configurations (this one included) copy the submitted script
  into a per-job spool directory before running it, which breaks any trick
  based on the script's own file location. `$SLURM_SUBMIT_DIR` is set
  correctly by SLURM regardless -- just make sure you `sbatch` from the
  same directory as your input files.
