#!/bin/bash
# Launch the Nextflow pipeline on CHPC SLURM. Run from the repo root on a login node.
#   ./run_nf.sh                 # real run (resumes if a previous run exists)
#   ./run_nf.sh -preview        # build the DAG, execute nothing
#   ./run_nf.sh --input other.csv
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

# Everything data-side is relative to --data_dir, so moving the storage (or
# running this on a laptop) is a one-flag change and samplesheet.csv stays put.
# --input is the exception: the sheet lives with the code, not with the data.
#
# References resolve against --db_dir, which is unset here and so falls back to
# --data_dir. When the databases move somewhere shared, add:
#     --db_dir /path/to/shared/databases \
# and drop the 'databases/' prefix from the two --sortmerna_* flags below.
#
# --sortmerna_idx reuses the index already built on this filesystem, skipping
# the ~16 min SORTMERNA_INDEX step. Drop the flag to rebuild it from scratch.
nextflow run "$HOMEDIR/main.nf" \
    -profile slurm \
    -work-dir "$SCRATCH/nf_work" \
    --data_dir "$SCRATCH" \
    --input "$HOMEDIR/samplesheet.csv" \
    --outdir "$SCRATCH/qc" \
    --sortmerna_ref databases/smr_v4.3_default_db.fasta \
    --sortmerna_idx databases/idx \
    -resume "$@"
