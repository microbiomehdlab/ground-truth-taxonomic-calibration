#!/usr/bin/env bash
# Prepare with the documented CLI first; this script submits only sealed plans.
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${ANALYSIS_SIF:?}"
concurrency="${REFERENCE_CONCURRENCY:-30}"
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid concurrency' >&2; exit 1; }
cd "$PROJECT"
tasks="$(python3 - "$REFERENCE_ROOT" "$ANALYSIS_SIF" <<'PY'
import csv,json,hashlib,sys
from pathlib import Path
root=Path(sys.argv[1]); state=json.loads((root/'identity.json').read_text())
if Path(sys.argv[2]).resolve()!=Path(state['image_path']): raise SystemExit('Image path differs')
for line in (root/'SHA256SUMS').read_text().splitlines():
    sha,name=line.split(None,1); p=(root/name).resolve()
    if root.resolve() not in p.parents or hashlib.sha256(p.read_bytes()).hexdigest()!=sha:
        raise SystemExit('Unsafe/changed plan member: '+name)
with (root/'tasks.tsv').open() as h: rows=list(csv.DictReader(h,delimiter='\t'))
if [int(r['index']) for r in rows]!=list(range(len(rows))) or not rows: raise SystemExit('Invalid tasks')
print(len(rows))
PY
)"
test ! -e "$REFERENCE_ROOT/jobs.tsv" || { echo 'Already submitted; inspect jobs.tsv before retrying' >&2; exit 1; }
job="$(sbatch --parsable --array="0-$((tasks-1))%$concurrency" --export=ALL analysis_v2/run_reference_response.sbatch)"
printf 'analysis\t%s\n' "$job" > "$REFERENCE_ROOT/jobs.tsv"
report="$(sbatch --parsable --dependency="afterok:${job%%;*}" --export=ALL analysis_v2/collect_reference_response.sbatch)"
printf 'report\t%s\n' "$report" >> "$REFERENCE_ROOT/jobs.tsv"
printf 'Analysis: %s\nReport: %s\nRoot: %s\n' "$job" "$report" "$REFERENCE_ROOT"
