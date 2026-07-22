#!/bin/bash
#SBATCH --job-name=brbseq_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
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
PIPELINE=EDIT_ME/brbseq-pipeline         # EDIT_ME: your actual GitHub owner/repo, e.g. cregan727/brbseq-pipeline
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
    --whitelist V5A_barcodes.txt \
    --outdir results \
    -resume
