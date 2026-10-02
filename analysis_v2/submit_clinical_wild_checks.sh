#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA_PILOT_RESULTS:?}" "${WILD_ROOT:?Fresh output required}"
cd "$PROJECT"
test ! -e "$WILD_ROOT"
test -s "$ANALYSIS_SIF"
test -d "$DA_PILOT_RESULTS"
[[ "${WILD_CONCURRENCY:-30}" =~ ^[1-9][0-9]*$ ]]
mkdir -p "$WILD_ROOT/source"
git archive HEAD | tar -x -C "$WILD_ROOT/source"
git rev-parse HEAD > "$WILD_ROOT/source_commit.txt"
export WILD_ROOT ANALYSIS_SIF DA_PILOT_RESULTS
prep="$(sbatch --parsable --job-name=wild_prepare --mem=4G --time=00:15:00 --output='wild_prepare.%j.out' --error='wild_prepare.%j.err' --export=ALL \
 --wrap='set -e; cd "$WILD_ROOT/source"; sha256sum "$ANALYSIS_SIF" > "$WILD_ROOT/image.sha256"; python3 -B analysis_v2/scripts/clinical_wild_checks.py plan --pilot "$DA_PILOT_RESULTS" --out "$WILD_ROOT/plan"')"
printf 'prepare\t%s\n' "$prep" | tee "$WILD_ROOT/jobs.tsv"
job="$(sbatch --parsable --array="0-155%${WILD_CONCURRENCY:-30}" --dependency="afterany:$prep" --export=ALL "$WILD_ROOT/source/analysis_v2/run_clinical_wild_checks.sbatch")"
printf 'validation\t%s\n' "$job" | tee -a "$WILD_ROOT/jobs.tsv"
report="$(sbatch --parsable --job-name=wild_report --dependency="afterany:$job" --mem=4G --time=00:15:00 --output='wild_report.%j.out' --error='wild_report.%j.err' --export=ALL \
 --wrap='set -e; cd "$WILD_ROOT/source"; python3 -B analysis_v2/scripts/clinical_wild_checks.py collect --plan "$WILD_ROOT/plan" --results "$WILD_ROOT/results" --out "$WILD_ROOT/REPORT"')"
printf 'report\t%s\n' "$report" | tee -a "$WILD_ROOT/jobs.tsv"
printf 'root\t%s\n' "$WILD_ROOT"
