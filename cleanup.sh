#!/bin/bash
# Retire the Snakemake pipeline and normalise the storage layout.
#
#   ./cleanup.sh        # show what would change, touch nothing
#   ./cleanup.sh --go   # actually do it
#
# Safe to re-run: every step checks whether it has already happened.
set -euo pipefail

REPO="$(cd "$(dirname "$0")" && pwd)"
DATA=${PVP_STORAGE:-/scratch/general/vast/abudge/plant_virus_pipeline_storage}

GO=0
case "${1:-}" in
    --go) GO=1 ;;
    ""|-n|--dry-run) GO=0 ;;
    *) echo "usage: $0 [--go]" >&2; exit 2 ;;
esac

(( GO )) || echo "== DRY RUN == nothing will be changed; re-run with --go"
echo

# --- refuse to run against a live pipeline -----------------------------------
# This renames the directory a running pipeline publishes into, which would
# strand its output half-written.
if squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -q '^nf-'; then
    echo "REFUSING: nextflow SLURM jobs are still running:" >&2
    squeue -u "$USER" -h -o '  %j %T %M' | grep '^  nf-' >&2
    exit 1
fi
if pgrep -u "$USER" -f 'nextflow run' >/dev/null 2>&1; then
    echo "REFUSING: a 'nextflow run' process is still active." >&2
    exit 1
fi

# --- refuse to delete qc/ unless qc_nf/ is a complete run --------------------
# qc/ is the Snakemake output. It is only disposable because the Nextflow run
# reproduced it; if that run is missing or partial, keep qc/.
if [[ -d "$DATA/qc" && -d "$DATA/qc_nf" ]]; then
    if ! compgen -G "$DATA/qc_nf/sortmerna/*/*.clean_fwd.fq.gz" >/dev/null \
       || [[ ! -f "$DATA/qc_nf/multiqc/multiqc_report.html" ]]; then
        echo "REFUSING: qc_nf/ has no clean reads or no multiqc report --" >&2
        echo "          it does not look like a completed run. Keeping qc/." >&2
        exit 1
    fi
fi

act() {  # act <description> <command...>
    local desc=$1; shift
    echo "  $desc"
    # 'if' rather than '&& ... || true': the latter swallows a failing rm.
    if (( GO )); then "$@"; fi
}

echo "--- repo: $REPO"
[[ -d "$REPO/attic" ]] \
    && act "rm  attic/ (retired Snakemake code; tarball on scratch)" rm -rf "$REPO/attic" \
    || echo "  ok  attic/ already gone"

shopt -s nullglob
rotated=( "$REPO"/.nextflow.log.[0-9]* )
if (( ${#rotated[@]} )); then
    act "rm  ${#rotated[@]} rotated .nextflow.log.N files (current log kept)" rm -f "${rotated[@]}"
else
    echo "  ok  no rotated nextflow logs"
fi
shopt -u nullglob

echo
echo "--- scratch: $DATA"
qc_gone=0
for d in qc:"Snakemake results (reproduced byte-for-byte by qc_nf)" \
         apptainer:"Snakemake apptainer-prefix" \
         apptainer_cache:"Snakemake APPTAINER_CACHEDIR" \
         _smrtest:"sortmerna SIGILL scratch test"; do
    name=${d%%:*}; why=${d#*:}
    if [[ -d "$DATA/$name" ]]; then
        act "rm  $name/  ($(du -sh "$DATA/$name" 2>/dev/null | cut -f1)) -- $why" rm -rf "$DATA/$name"
        [[ $name == qc ]] && qc_gone=1
    else
        echo "  ok  $name/ already gone"
        [[ $name == qc ]] && qc_gone=1
    fi
done

if [[ -d "$DATA/qc_nf" ]] && (( qc_gone )); then
    act "mv  qc_nf/ -> qc/  (Nextflow output takes the canonical name)" mv "$DATA/qc_nf" "$DATA/qc"
elif [[ -d "$DATA/qc_nf" ]]; then
    echo "  !!  qc/ still present and not scheduled for removal -- not renaming"
else
    echo "  ok  qc_nf/ already renamed"
fi

# --- point run.sh at the canonical outdir ------------------------------------
if [[ -f "$REPO/run.sh" ]] && grep -q 'qc_nf' "$REPO/run.sh"; then
    act "ed  run.sh: --outdir qc_nf -> qc, drop the stale qc_nf comment" \
        python3 - "$REPO/run.sh" <<'PY'
import re, sys, pathlib
p = pathlib.Path(sys.argv[1]); s = p.read_text()
s = s.replace('--outdir "$SCRATCH/qc_nf"', '--outdir "$SCRATCH/qc"')
s = re.sub(r'#\n# --outdir is still qc_nf because.*?--outdir "\$SCRATCH/qc"\.\n',
           '', s, flags=re.S)
p.write_text(s)
PY
else
    echo "  ok  run.sh already points at qc/"
fi

echo
echo "--- preserved (not touched)"
for keep in raw:"input reads" \
            databases:"sortmerna reference + prebuilt index" \
            fastas:"uniref90 / NCBI viral -- staged for future processes" \
            apptainer_nf:"Nextflow container cache" \
            nf_work:"Nextflow work dir (reclaim with: nextflow clean -f -before <run>)"; do
    name=${keep%%:*}; why=${keep#*:}
    [[ -e "$DATA/$name" ]] && printf "  %-14s %-6s %s\n" "$name/" \
        "$(du -sh "$DATA/$name" 2>/dev/null | cut -f1)" "$why"
done
ls "$DATA"/snakemake_retired_*.tar.gz >/dev/null 2>&1 \
    && echo "  snakemake_retired_*.tar.gz   Snakemake source backup"

echo
(( GO )) && echo "done." || echo "dry run only -- re-run with --go to apply."
