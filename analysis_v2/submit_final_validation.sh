#!/bin/bash
# One entry point; queues all checks and an after-any collector.
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA2_PAIRED_RESULTS:?}" "${DA_INVENTORY_ROOT:?}"
cd "$PROJECT"
export FINAL_VALIDATION_ROOT="${FINAL_VALIDATION_ROOT:-$PROJECT/work/final_da_validation_$(date -u +%Y%m%dT%H%M%SZ)}"
test ! -e "$FINAL_VALIDATION_ROOT"
test -s "$ANALYSIS_SIF"
test -d "$DA2_PAIRED_RESULTS"
test -d "$DA_INVENTORY_ROOT"
# Enforce actual backend comparison in the pinned container, not a local skip.
apptainer exec --cleanenv --env REQUIRE_VALIDATION_BACKEND=1 \
  --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
  Rscript analysis_v2/tests/test_final_validation.R
mkdir -p "$FINAL_VALIDATION_ROOT"
git rev-parse HEAD > "$FINAL_VALIDATION_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" analysis_v2/scripts/run_final_validation_task.py \
  analysis_v2/scripts/run_final_validation.R analysis_v2/scripts/collect_final_validation.py \
  analysis_v2/lib/unpaired_null.R analysis_v2/lib/paired_null.R \
  analysis_v2/lib/paired_difference_context.R analysis_v2/lib/maaslin_contract.R \
  analysis_v2/lib/maaslin_context.R analysis_v2/run_final_validation.sbatch \
  analysis_v2/scripts/prepare_da_pilot.py analysis_v2/scripts/audit_bracken_denominators.py \
  analysis_v2/scripts/build_biomarker_abundance_input.py \
  analysis_v2/collect_final_validation.sbatch > "$FINAL_VALIDATION_ROOT/submission.sha256"
validation_job="$(sbatch --parsable --array=0-39%${VALIDATION_CONCURRENCY:-30} \
  --export=ALL analysis_v2/run_final_validation.sbatch)"
report_job="$(sbatch --parsable --dependency="afterany:${validation_job%%;*}" \
  --export=ALL analysis_v2/collect_final_validation.sbatch)"
printf 'validation_job\t%s\nreport_job\t%s\nroot\t%s\n' \
  "$validation_job" "$report_job" "$FINAL_VALIDATION_ROOT" | tee "$FINAL_VALIDATION_ROOT/jobs.tsv"
