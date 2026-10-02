#!/bin/bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${FINAL_VALIDATION_ROOT:?}" "${DA_PILOT_RESULTS:?}" "${DA_PILOT_PLAN:?}"
cd "$PROJECT"
export CANDIDATE_ROOT="${CANDIDATE_ROOT:-$PROJECT/work/da3_exact_candidate_$(date -u +%Y%m%dT%H%M%SZ)}"
test ! -e "$CANDIDATE_ROOT"
export CANDIDATE_SOURCE="$CANDIDATE_ROOT/source"
mkdir -p "$CANDIDATE_SOURCE"
# Immutable tracked source snapshot lets other development/pulls proceed safely.
git archive HEAD | tar -x -C "$CANDIDATE_SOURCE"
git rev-parse HEAD > "$CANDIDATE_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$CANDIDATE_ROOT/analysis_image.sha256"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$CANDIDATE_SOURCE" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_exact_da3_candidate.R
cd "$CANDIDATE_SOURCE"
unset CANDIDATE_COLLECT
job="$(sbatch --parsable --array=0-51%${VALIDATION_CONCURRENCY:-30} --export=ALL analysis_v2/run_da3_candidate.sbatch)"
export CANDIDATE_COLLECT=1
report="$(sbatch --parsable --dependency="afterany:${job%%;*}" --export=ALL \
  --job-name=da3_candidate_report --output=da3_candidate_report.%j.out --error=da3_candidate_report.%j.err \
  analysis_v2/run_da3_candidate.sbatch)"
printf 'candidate_job\t%s\nreport_job\t%s\nroot\t%s\n' "$job" "$report" "$CANDIDATE_ROOT" | tee "$CANDIDATE_ROOT/jobs.tsv"
