#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA_PILOT_RESULTS:?}" "${DA_INVENTORY_ROOT:?}" "${CHECK_ROOT:?Use fresh output}"
cd "$PROJECT"
test ! -e "$CHECK_ROOT"
test -s "$ANALYSIS_SIF"
mkdir -p "$CHECK_ROOT/source"
# Use committed immutable code, not the live working tree or unrelated edits.
git archive HEAD | tar -x -C "$CHECK_ROOT/source"
git rev-parse HEAD > "$CHECK_ROOT/source_commit.txt"
export CHECK_ROOT ANALYSIS_SIF DA_PILOT_RESULTS DA_INVENTORY_ROOT
prep="$(sbatch --parsable --job-name=preprod_prepare --mem=8G --time=00:30:00 \
    --output='preproduction_prepare.%j.out' --error='preproduction_prepare.%j.err' --export=ALL \
    --wrap='set -e; cd "$CHECK_ROOT/source"; sha256sum "$ANALYSIS_SIF" > "$CHECK_ROOT/analysis_image.sha256"; python3 -B analysis_v2/scripts/remaining_validation.py plan --pilot "$DA_PILOT_RESULTS" --out "$CHECK_ROOT/validation_plan"; python3 -B analysis_v2/scripts/plan_da3_batches.py plan --inventory "$DA_INVENTORY_ROOT" --out "$CHECK_ROOT/draft_da3_plan"')"
export CHECK_MODE=validation
validation="$(sbatch --parsable --array="0-71%${CHECK_CONCURRENCY:-15}" --dependency="afterany:$prep" --export=ALL "$CHECK_ROOT/source/analysis_v2/run_preproduction_task.sbatch")"
export CHECK_MODE=batch_canary
canary="$(sbatch --parsable --array="0-47%${CHECK_CONCURRENCY:-15}" --dependency="afterany:$prep" --export=ALL "$CHECK_ROOT/source/analysis_v2/run_preproduction_task.sbatch")"
collector="$(sbatch --parsable --job-name=preprod_report --mem=4G --time=00:15:00 \
    --dependency="afterany:$validation:$canary" --output='preproduction_report.%j.out' --error='preproduction_report.%j.err' --export=ALL \
    --wrap='cd "$CHECK_ROOT/source"; python3 -B analysis_v2/scripts/collect_preproduction.py "$CHECK_ROOT"')"
printf 'prepare\t%s\nvalidation\t%s\ncompact_canary\t%s\nreport\t%s\nroot\t%s\n' "$prep" "$validation" "$canary" "$collector" "$CHECK_ROOT"
