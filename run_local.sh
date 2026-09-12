#!/bin/bash
# Run the pipeline locally through apptainer for fast development iteration,
# against a small subset of reads. Uses the same storage layout as the HPC run.
#   PVP_STORAGE=~/plant_virus_data ./run_local.sh
#
# Data layout under $PVP_STORAGE (mirrors CHPC):
#   raw/           reads — a subset is fine
#   databases/     fixed database layout documented in README.md
#   apptainer_nf/  (optional) copy the HPC cache dir here to skip pulling;
#                  otherwise apptainer pulls from quay.io on first use.
set -euo pipefail

HOMEDIR=$(cd "$(dirname "$0")" && pwd)
PVP_STORAGE=${PVP_STORAGE:-$HOME/plant_virus_data}

export NXF_APPTAINER_CACHEDIR="$PVP_STORAGE/apptainer_nf"
mkdir -p "$NXF_APPTAINER_CACHEDIR"

nextflow run "$HOMEDIR/main.nf" \
    -profile local \
    -work-dir "$PVP_STORAGE/nf_work" \
    --data_dir "$PVP_STORAGE" \
    --input "$HOMEDIR/samplesheet.tsv" \
    --outdir "$PVP_STORAGE/qc" \
    --db databases \
    "$@"
