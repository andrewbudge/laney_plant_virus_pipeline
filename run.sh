#!/bin/bash
# Launch the plant-virus QC pipeline on CHPC SLURM. Run from the repo root on a
# login node.
#   ./run.sh               # fresh real run
#   ./run.sh -preview      # build the DAG, execute nothing
set -euo pipefail

module load nextflow/25.04
module load apptainer/1.4.0

# CHPC storage root. Override on another machine:
#   PVP_STORAGE=~/plant_virus_data ./run.sh
SCRATCH=${PVP_STORAGE:-/scratch/general/vast/abudge/plant_virus_pipeline_storage}
HOMEDIR=$(cd "$(dirname "$0")" && pwd)

# Keep apptainer's pull cache off $HOME — it has a quota, scratch does not.
export NXF_APPTAINER_CACHEDIR="$SCRATCH/apptainer_nf"
mkdir -p "$NXF_APPTAINER_CACHEDIR"

# FASTQ paths are relative to --data_dir. All references use the fixed --db
# layout documented in README.md.
nextflow run "$HOMEDIR/main.nf" \
    -profile slurm \
    -work-dir "$SCRATCH/nf_work" \
    --data_dir "$SCRATCH" \
    --input "$HOMEDIR/samplesheet.tsv" \
    --outdir "$SCRATCH/qc" \
    --db databases \
    "$@"
