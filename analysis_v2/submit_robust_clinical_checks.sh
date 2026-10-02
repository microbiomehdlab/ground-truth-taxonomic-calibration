#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA_PILOT_RESULTS:?}" "${CHECK_ROOT:?Fresh output required}"
cd "$PROJECT";test ! -e "$CHECK_ROOT";test -s "$ANALYSIS_SIF"
mkdir -p "$CHECK_ROOT/source"
git archive HEAD | tar -x -C "$CHECK_ROOT/source"
git rev-parse HEAD > "$CHECK_ROOT/source_commit.txt"
export CHECK_ROOT ANALYSIS_SIF DA_PILOT_RESULTS CHECK_MODE=validation
prep="$(sbatch --parsable --job-name=hc3_prepare --mem=4G --time=00:15:00 \
 --output='hc3_prepare.%j.out' --error='hc3_prepare.%j.err' --export=ALL \
 --wrap='set -e; cd "$CHECK_ROOT/source"; sha256sum "$ANALYSIS_SIF" > "$CHECK_ROOT/image.sha256"; python3 -B analysis_v2/scripts/remaining_validation.py plan --clinical-only --pilot "$DA_PILOT_RESULTS" --out "$CHECK_ROOT/validation_plan"')"
job="$(sbatch --parsable --array="0-47%${CHECK_CONCURRENCY:-30}" --dependency="afterany:$prep" --export=ALL "$CHECK_ROOT/source/analysis_v2/run_preproduction_task.sbatch")"
report="$(sbatch --parsable --job-name=hc3_report --dependency="afterany:$job" --mem=4G --time=00:15:00 \
 --output='hc3_report.%j.out' --error='hc3_report.%j.err' --export=ALL \
 --wrap='cd "$CHECK_ROOT/source"; python3 -B analysis_v2/scripts/remaining_validation.py collect --plan "$CHECK_ROOT/validation_plan" --results "$CHECK_ROOT/validation_results" --out "$CHECK_ROOT/REPORT"')"
printf 'prepare\t%s\nvalidation\t%s\nreport\t%s\nroot\t%s\n' "$prep" "$job" "$report" "$CHECK_ROOT"
