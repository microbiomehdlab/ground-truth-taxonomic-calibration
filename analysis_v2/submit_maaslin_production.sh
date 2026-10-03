#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${MAASLIN_RUN_ROOT:?}"
cd "$PROJECT"
mode="${1:-new}"
[[ "$mode" == new || "$mode" == --resume ]]
test -s "$ANALYSIS_SIF"
export MAASLIN_CONCURRENCY="${MAASLIN_CONCURRENCY:-30}"
if ! [[ "$MAASLIN_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] || ((MAASLIN_CONCURRENCY<6 || MAASLIN_CONCURRENCY>192)); then
    echo 'MAASLIN_CONCURRENCY must be an integer from 6 to 192' >&2
    exit 1
fi
if [[ "$mode" == new ]]; then
    : "${DA3_BATCH_PLAN:?}"
    python3 - "$MAASLIN_RUN_ROOT" "$DA3_BATCH_PLAN" "$PROJECT" <<'PY'
import sys
from pathlib import Path
out,plan,project=[Path(x).resolve() for x in sys.argv[1:]]
if out.exists() or out==plan or out in plan.parents or plan in out.parents or out==project or out in project.parents:
    raise SystemExit('Use a fresh output root, disjoint from the plan and not an ancestor of the project')
PY
    mkdir "$MAASLIN_RUN_ROOT"
    mkdir "$MAASLIN_RUN_ROOT/source"
    git archive HEAD | tar -x -C "$MAASLIN_RUN_ROOT/source"
    git rev-parse HEAD > "$MAASLIN_RUN_ROOT/source_commit.txt"
    export MAASLIN_RESUME=0
else
    test -s "$MAASLIN_RUN_ROOT/identity.json"
    export MAASLIN_RESUME=1
fi
exec 9>"$MAASLIN_RUN_ROOT/submission.lock"
flock -n 9 || { echo 'Another submission is in progress' >&2; exit 1; }
if [[ "$mode" == --resume ]]; then
    # Do not resubmit while any job previously recorded for this root is active.
    active="$(squeue -h -u "$USER" -o '%i')"
    python3 - "$MAASLIN_RUN_ROOT/jobs.tsv" "$active" <<'PY'
import sys
from pathlib import Path
path=Path(sys.argv[1])
old={line.split('\t')[1].split(';')[0] for line in path.read_text().splitlines()} if path.exists() else set()
active={line.strip().split('_')[0] for line in sys.argv[2].splitlines()}
if old & active: raise SystemExit('Recorded jobs still active; do not resume yet: '+str(sorted(old & active)))
PY
fi
export MAASLIN_RUN_ROOT ANALYSIS_SIF
script="$MAASLIN_RUN_ROOT/source/analysis_v2/run_maaslin_production.sbatch"
test -s "$script"
export MAASLIN_REPORT_ROOT="$MAASLIN_RUN_ROOT/REPORT_$(date -u +%Y%m%dT%H%M%SZ)"
test ! -e "$MAASLIN_REPORT_ROOT"
prep="$(sbatch --parsable --export=ALL --job-name=da3_maaslin_prepare \
    --output='da3_maaslin_prepare.%j.out' --error='da3_maaslin_prepare.%j.err' "$script" prepare)"
printf 'prepare\t%s\n' "$prep" | tee -a "$MAASLIN_RUN_ROOT/jobs.tsv"
dependencies=""
for chunk in 0 1 2 3 4 5; do
    export MAASLIN_BATCH_OFFSET=$((chunk*800))
    throttle=$((MAASLIN_CONCURRENCY/6 + (chunk < MAASLIN_CONCURRENCY%6)))
    job="$(sbatch --parsable --export=ALL --array="0-799%$throttle" \
        --dependency="afterok:${prep%%;*}" "$script" run)"
    printf 'array_%s\t%s\n' "$MAASLIN_BATCH_OFFSET" "$job" | tee -a "$MAASLIN_RUN_ROOT/jobs.tsv"
    dependencies="${dependencies:+$dependencies:}${job%%;*}"
done
report="$(sbatch --parsable --export=ALL --job-name=da3_maaslin_collect --time=06:00:00 \
    --dependency="afterany:$dependencies" --output='da3_maaslin_collect.%j.out' \
    --error='da3_maaslin_collect.%j.err' "$script" collect)"
printf 'collect\t%s\n' "$report" | tee -a "$MAASLIN_RUN_ROOT/jobs.tsv"
printf '%s\n' "$MAASLIN_REPORT_ROOT" >> "$MAASLIN_RUN_ROOT/report_paths.txt"
printf 'root\t%s\nreport\t%s\nconcurrent_workers\t%s\n' "$MAASLIN_RUN_ROOT" "$MAASLIN_REPORT_ROOT" "$MAASLIN_CONCURRENCY"
