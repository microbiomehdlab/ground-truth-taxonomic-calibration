#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}" "${MAASLIN_TIMING_ROOT:?Fresh root required}"
cd "$PROJECT"
test -s "$ANALYSIS_SIF"
test -s "$DA3_BATCH_PLAN/SHA256SUMS"
export MAASLIN_TIMING_CONCURRENCY="${MAASLIN_TIMING_CONCURRENCY:-12}"
[[ "$MAASLIN_TIMING_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] && ((MAASLIN_TIMING_CONCURRENCY<=48))
python3 - "$MAASLIN_TIMING_ROOT" "$DA3_BATCH_PLAN" "$PROJECT" <<'PY'
from pathlib import Path
import sys
out, plan, project = [Path(p).resolve() for p in sys.argv[1:]]
if out.exists() or out==plan or plan in out.parents or out in plan.parents or out==project or out in project.parents:
    raise SystemExit('Choose a fresh output root that cannot overlap the input plan or contain the project')
PY
mkdir -p "$MAASLIN_TIMING_ROOT/source"
git archive HEAD | tar -x -C "$MAASLIN_TIMING_ROOT/source"
git rev-parse HEAD > "$MAASLIN_TIMING_ROOT/source_commit.txt"
test -s "$MAASLIN_TIMING_ROOT/source/analysis_v2/scripts/maaslin_timing.py"
export MAASLIN_TIMING_ROOT DA3_BATCH_PLAN ANALYSIS_SIF
script="$MAASLIN_TIMING_ROOT/source/analysis_v2/run_maaslin_timing.sbatch"
prep="$(sbatch --parsable --export=ALL --job-name=maaslin_prepare \
    --output='maaslin_prepare.%j.out' --error='maaslin_prepare.%j.err' "$script" prepare)"
printf 'prepare\t%s\n' "$prep" | tee "$MAASLIN_TIMING_ROOT/jobs.tsv"
job="$(sbatch --parsable --export=ALL --array="0-47%$MAASLIN_TIMING_CONCURRENCY" \
    --dependency="afterany:${prep%%;*}" "$script" run)"
printf 'timing\t%s\n' "$job" | tee -a "$MAASLIN_TIMING_ROOT/jobs.tsv"
report="$(sbatch --parsable --export=ALL --job-name=maaslin_report \
    --dependency="afterany:${job%%;*}" --output='maaslin_report.%j.out' \
    --error='maaslin_report.%j.err' "$script" collect)"
printf 'report\t%s\n' "$report" | tee -a "$MAASLIN_TIMING_ROOT/jobs.tsv"
printf 'root\t%s\n' "$MAASLIN_TIMING_ROOT"
