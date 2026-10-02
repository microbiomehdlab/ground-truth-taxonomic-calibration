#!/bin/bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${FINAL_VALIDATION_ROOT:?}"
cd "$PROJECT"
export ROBUSTNESS_ROOT="${ROBUSTNESS_ROOT:-$PROJECT/work/da_robustness_$(date -u +%Y%m%dT%H%M%SZ)}"
test ! -e "$ROBUSTNESS_ROOT"
export ROBUSTNESS_SOURCE="$ROBUSTNESS_ROOT/source"
mkdir -p "$ROBUSTNESS_SOURCE"
git archive HEAD | tar -x -C "$ROBUSTNESS_SOURCE"
git rev-parse HEAD > "$ROBUSTNESS_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$ROBUSTNESS_ROOT/analysis_image.sha256"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$ROBUSTNESS_SOURCE" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_da_robustness.R
cd "$ROBUSTNESS_SOURCE"
unset ROBUSTNESS_COLLECT
job="$(sbatch --parsable --array=0-387%${VALIDATION_CONCURRENCY:-30} --export=ALL analysis_v2/run_da_robustness.sbatch)"
export ROBUSTNESS_COLLECT=1
report="$(sbatch --parsable --dependency="afterany:${job%%;*}" --export=ALL --job-name=da_robustness_report \
  --output=da_robustness_report.%j.out --error=da_robustness_report.%j.err analysis_v2/run_da_robustness.sbatch)"
printf 'robustness_job\t%s\nreport_job\t%s\nroot\t%s\n' "$job" "$report" "$ROBUSTNESS_ROOT" | tee "$ROBUSTNESS_ROOT/jobs.tsv"
