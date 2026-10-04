#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA_INVENTORY_ROOT:?}"
: "${BRACKEN_RECOVERY_ROOT:?}" "${METAPHLAN_RECOVERY_ROOT:?}" "${DA2_ROOT:?}"
export DA2_CONCURRENCY="${DA2_CONCURRENCY:-6}"
[[ "$DA2_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] && ((DA2_CONCURRENCY<=18))
cd "$PROJECT"
python3 - "$DA2_ROOT" "$PROJECT" "$DA_INVENTORY_ROOT" "$BRACKEN_RECOVERY_ROOT" "$METAPHLAN_RECOVERY_ROOT" <<'PY'
from pathlib import Path
import sys
out,project,*inputs=[Path(p).resolve() for p in sys.argv[1:]]
assert not out.exists() and out!=project and out not in project.parents, 'Fresh output required'
assert all(out!=p and out not in p.parents and p not in out.parents for p in inputs), 'Overlapping output'
PY
test -s "$ANALYSIS_SIF"
mkdir -p "$DA2_ROOT/source"
git archive HEAD | tar -x -C "$DA2_ROOT/source"
git rev-parse HEAD > "$DA2_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$DA2_ROOT/image.sha256"
(cd "$DA2_ROOT/source" && find . -type f ! -name SOURCE_SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SOURCE_SHA256SUMS)
export PROJECT ANALYSIS_SIF DA2_ROOT DA_INVENTORY_ROOT BRACKEN_RECOVERY_ROOT METAPHLAN_RECOVERY_ROOT
script="$DA2_ROOT/source/analysis_v2/run_da2_dose_response.sbatch"
prep="$(sbatch --parsable --export=ALL "$script" prepare)"
printf 'prepare\t%s\n' "$prep" | tee "$DA2_ROOT/jobs.tsv"
job="$(sbatch --parsable --export=ALL --dependency="afterany:${prep%%;*}" --array="0-17%$DA2_CONCURRENCY" "$script" run)"
printf 'analysis\t%s\n' "$job" | tee -a "$DA2_ROOT/jobs.tsv"
report="$(sbatch --parsable --export=ALL --dependency="afterany:${job%%;*}" \
    --output='da2_report.%j.out' --error='da2_report.%j.err' "$script" collect)"
printf 'report\t%s\n' "$report" | tee -a "$DA2_ROOT/jobs.tsv"
printf 'root\t%s\n' "$DA2_ROOT"
