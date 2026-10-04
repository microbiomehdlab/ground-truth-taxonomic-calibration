#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA1_REPORT:?}" "${DA_INVENTORY_ROOT:?}" "${CLINICAL_SPIKE_ROOT:?}"
export CLINICAL_SPIKE_CONCURRENCY="${CLINICAL_SPIKE_CONCURRENCY:-30}"
[[ "$CLINICAL_SPIKE_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] && ((CLINICAL_SPIKE_CONCURRENCY<=42))
cd "$PROJECT"
python3 - "$CLINICAL_SPIKE_ROOT" "$PROJECT" "$DA1_REPORT" "$DA_INVENTORY_ROOT" <<'PY'
from pathlib import Path
import sys
out, project, *inputs = [Path(p).resolve() for p in sys.argv[1:]]
assert not out.exists() and out != project and out not in project.parents, 'Fresh output required'
assert all(out != p and out not in p.parents and p not in out.parents for p in inputs), 'Output overlaps input'
PY
git cat-file -e HEAD:analysis_v2/scripts/clinical_spike_association.py
test -s "$ANALYSIS_SIF"
mkdir -p "$CLINICAL_SPIKE_ROOT/source"
git archive HEAD | tar -x -C "$CLINICAL_SPIKE_ROOT/source"
git rev-parse HEAD > "$CLINICAL_SPIKE_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$CLINICAL_SPIKE_ROOT/image.sha256"
(cd "$CLINICAL_SPIKE_ROOT/source" && find . -type f ! -name SOURCE_SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SOURCE_SHA256SUMS)
script="$CLINICAL_SPIKE_ROOT/source/analysis_v2/run_clinical_spike.sbatch"
prep="$(sbatch --parsable --export=ALL "$script" prepare)"
printf 'prepare\t%s\n' "$prep" | tee "$CLINICAL_SPIKE_ROOT/jobs.tsv"
job="$(sbatch --parsable --export=ALL --dependency="afterok:${prep%%;*}" --array="0-41%$CLINICAL_SPIKE_CONCURRENCY" "$script" run)"
printf 'analysis\t%s\n' "$job" | tee -a "$CLINICAL_SPIKE_ROOT/jobs.tsv"
report="$(sbatch --parsable --export=ALL --dependency="afterany:${job%%;*}" --output='clinical_spike_report.%j.out' --error='clinical_spike_report.%j.err' "$script" collect)"
printf 'report\t%s\n' "$report" | tee -a "$CLINICAL_SPIKE_ROOT/jobs.tsv"
printf 'root\t%s\n' "$CLINICAL_SPIKE_ROOT"
