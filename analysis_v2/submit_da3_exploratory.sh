#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}" "${DA3_RUN_ROOT:?Fresh output required}"
cd "$PROJECT"
test ! -e "$DA3_RUN_ROOT"
test -s "$ANALYSIS_SIF"
# Read-only input verification before creating the run or submitting any jobs.
python3 -B analysis_v2/scripts/plan_da3_batches.py check --plan "$DA3_BATCH_PLAN"
mkdir -p "$DA3_RUN_ROOT/source"
git archive HEAD | tar -x -C "$DA3_RUN_ROOT/source"
git rev-parse HEAD > "$DA3_RUN_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$DA3_RUN_ROOT/analysis_image.sha256"
read -r DA3_IMAGE_SHA256 image_path < "$DA3_RUN_ROOT/analysis_image.sha256"
export DA3_IMAGE_SHA256
export DA3_BATCH_SOURCE="$DA3_RUN_ROOT/source"
export DA3_BATCH_RESULTS="$DA3_RUN_ROOT/results"
export DA3_BATCH_PLAN ANALYSIS_SIF
jobs=()
for offset in 0 800 1600 2400 3200 4000; do
    export DA3_BATCH_OFFSET="$offset"
    # Six arrays × five concurrent workers = at most thirty workers overall.
    job="$(sbatch --parsable --array=0-799%5 --export=ALL "$DA3_BATCH_SOURCE/analysis_v2/run_da3_batches.sbatch")"
    jobs+=("$job")
    printf '%s\t%s\n' "$offset" "$job" >> "$DA3_RUN_ROOT/submitted_arrays.tsv"
done
dependency="$(IFS=:; echo "${jobs[*]}")"
report="$(sbatch --parsable --job-name=da3_collect --dependency="afterany:$dependency" \
    --mem=8G --time=00:30:00 --output='da3_collect.%j.out' --error='da3_collect.%j.err' --export=ALL \
    --wrap='cd "$DA3_BATCH_SOURCE"; python3 -B analysis_v2/scripts/collect_da3_batches.py --plan "$DA3_BATCH_PLAN" --results "$DA3_BATCH_RESULTS" --out "$DA3_BATCH_RESULTS/../REPORT"')"
printf 'arrays\t%s\ncollector\t%s\nroot\t%s\n' "${jobs[*]}" "$report" "$DA3_RUN_ROOT"
