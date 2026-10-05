#!/usr/bin/env bash
#SBATCH --job-name=reference_prepare
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=reference_prepare.%j.out
#SBATCH --error=reference_prepare.%j.err
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}"
: "${DA_INVENTORY_ROOT:?}" "${GEFF_ROOT:?}" "${DA3_OBSERVED_REPORT:?}"
allocations="${REFERENCE_ALLOCATIONS:-20}"
[[ "$allocations" == 20 || "$allocations" == 100 ]] || { echo 'Invalid REFERENCE_ALLOCATIONS' >&2; exit 1; }
cd "$PROJECT"
runtime=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --bind /mnt/beegfs:/mnt/beegfs
    --pwd "$PROJECT" --env TZ=UTC,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 "$ANALYSIS_SIF")
"${runtime[@]}" python3 -B analysis_v2/tests/test_reference_response_pilot.py
"${runtime[@]}" Rscript analysis_v2/tests/test_reference_response.R "$PROJECT" --actual
"${runtime[@]}" python3 -B analysis_v2/scripts/reference_response_pilot.py prepare \
    --plan "$DA3_BATCH_PLAN" --inventory "$DA_INVENTORY_ROOT" --geff-root "$GEFF_ROOT" \
    --observed-report "$DA3_OBSERVED_REPORT" --repo "$PROJECT" --image "$ANALYSIS_SIF" \
    --out "$REFERENCE_ROOT" --allocations "$allocations" --arms U --batch-size 10
python3 - "$REFERENCE_ROOT" "$allocations" <<'PY'
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1]); expected=36*int(sys.argv[2])
with (root/'tasks.tsv').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
assert len(rows)==expected//10
assert [int(r['index']) for r in rows]==list(range(len(rows)))
assert all(int(r['contexts'])==10 for r in rows)
assert json.loads((root/'identity.json').read_text())['contexts']==expected
print('[PASS] Actual pinned backend and immutable reference plan:',expected,'contexts,',len(rows),'tasks')
PY
