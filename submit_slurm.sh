#!/bin/bash
#SBATCH --job-name=brbseq_pipeline
# Sized to cover STARSOLO's per-task request (cpus 32, memory 120 GB in
# nextflow.config) plus headroom for other steps running concurrently --
# NOT just "enough to launch Nextflow itself". executor = 'local' runs
# every process inside THIS allocation (no sub-job submission), so
# undersizing this doesn't fail cleanly -- STARSOLO gets silently
# OOM-killed by the cgroup partway through a 48-hour job instead. Keep
# this in sync with the `executor { cpus / memory }` block in
# nextflow.config; bump both together if you want more than one STARSOLO
# task (one per library x genome pair) to run at once.
#SBATCH --cpus-per-task=32
#SBATCH --mem=140G
#SBATCH --time=48:00:00
#SBATCH --output=logs/brbseq_%j.out
#SBATCH --error=logs/brbseq_%j.err

set -euo pipefail

# Run `sbatch` from this same directory (the one containing main.nf,
# nextflow.config, and your input files) -- $SLURM_SUBMIT_DIR is set
# correctly by SLURM regardless of submission directory, unlike path tricks
# based on this script's own file location, which break on clusters (like
# the one this was developed against) that copy the submitted script into a
# per-job spool directory before running it.
SCRIPT_DIR="$SLURM_SUBMIT_DIR"
mkdir -p "$SCRIPT_DIR/logs"
cd "$SCRIPT_DIR"

### ---- EDIT THESE ----------------------------------------------------------
NEXTFLOW_BIN=nextflow                    # EDIT_ME: path to your nextflow binary/module, e.g. a conda env
PIPELINE=cregan727/nf-rnator              # your actual GitHub owner/repo
PIPELINE_REVISION=main                   # branch, tag, or commit -- pin this once the pipeline is stable
module load singularity                  # EDIT_ME: check `module avail singularity apptainer` and match nextflow.config
CONTAINER_PROFILE=singularity            # change if your cluster uses apptainer/docker instead
### ---------------------------------------------------------------------------

export NXF_SINGULARITY_CACHEDIR="$SCRIPT_DIR/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

"$NEXTFLOW_BIN" run "$PIPELINE" -r "$PIPELINE_REVISION" \
    -profile "$CONTAINER_PROFILE" \
    --input samplesheet.csv \
    --genomes genomes.csv \
    --whitelist V5D_barcodes.txt \
    --outdir results \
    -resume
