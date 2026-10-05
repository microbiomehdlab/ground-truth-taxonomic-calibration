#!/usr/bin/env bash
# One engineering submission; preparation gates both the array and collector.
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}"
: "${DA_INVENTORY_ROOT:?}" "${GEFF_ROOT:?}" "${DA3_OBSERVED_REPORT:?}"
concurrency="${REFERENCE_CONCURRENCY:-30}"
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid concurrency' >&2; exit 1; }
test ! -e "$REFERENCE_ROOT" && test ! -e "$REFERENCE_ROOT.jobs.tsv" || {
    echo 'Use a fresh REFERENCE_ROOT; no automatic resubmission' >&2; exit 1;
}
cd "$PROJECT"
for path in "$DA3_BATCH_PLAN/SHA256SUMS" "$DA_INVENTORY_ROOT/yachida/SHA256SUMS" \
    "$GEFF_ROOT/SUCCESS" "$DA3_OBSERVED_REPORT/status.json" "$ANALYSIS_SIF"; do
    test -s "$path" || { printf 'Missing input: %s\n' "$path" >&2; exit 1; }
done
prepare="$(sbatch --parsable --export=ALL analysis_v2/prepare_reference_response.sbatch)"
printf 'prepare\t%s\n' "$prepare" > "$REFERENCE_ROOT.jobs.tsv"
job="$(sbatch --parsable --array="0-71%$concurrency" --dependency="afterok:${prepare%%;*}" \
    --export=ALL analysis_v2/run_reference_response.sbatch)"
printf 'analysis\t%s\n' "$job" >> "$REFERENCE_ROOT.jobs.tsv"
report="$(sbatch --parsable --dependency="afterok:${job%%;*}" --export=ALL analysis_v2/collect_reference_response.sbatch)"
printf 'report\t%s\n' "$report" >> "$REFERENCE_ROOT.jobs.tsv"
printf 'Preparation: %s\nAnalysis: %s\nReport: %s\nRoot: %s\n' "$prepare" "$job" "$report" "$REFERENCE_ROOT"
