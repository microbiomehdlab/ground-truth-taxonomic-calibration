#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${REFERENCE_SCALE_ROOT:?}" "${ANALYSIS_SIF:?}"
concurrency="${REFERENCE_CONCURRENCY:-30}"
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid concurrency' >&2; exit 1; }
test ! -e "$REFERENCE_SCALE_ROOT" && test ! -e "$REFERENCE_SCALE_ROOT.jobs.tsv" || {
    echo 'Use a fresh sensitivity root; never overwrite completed runs' >&2; exit 1;
}
cd "$PROJECT"
prepare="$(sbatch --parsable --export=ALL analysis_v2/prepare_reference_scaling.sbatch)"
printf 'prepare\t%s\n' "$prepare" > "$REFERENCE_SCALE_ROOT.jobs.tsv"
analysis="$(sbatch --parsable --export=ALL --array="0-71%$concurrency" \
    --dependency="afterok:${prepare%%;*}" analysis_v2/run_reference_scaling.sbatch)"
printf 'analysis\t%s\n' "$analysis" >> "$REFERENCE_SCALE_ROOT.jobs.tsv"
report="$(sbatch --parsable --export=ALL --dependency="afterok:${analysis%%;*}" analysis_v2/collect_reference_scaling.sbatch)"
printf 'report\t%s\n' "$report" >> "$REFERENCE_SCALE_ROOT.jobs.tsv"
printf 'Prepare: %s\nAnalysis: %s\nReport: %s\nOutput: %s\n' "$prepare" "$analysis" "$report" "$REFERENCE_SCALE_ROOT"
