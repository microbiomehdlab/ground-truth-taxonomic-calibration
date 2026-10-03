#!/usr/bin/env bash
#SBATCH --job-name=da3_maaslin
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=da3_maaslin.%A_%a.out
#SBATCH --error=da3_maaslin.%A_%a.err
set -euo pipefail
: "${MAASLIN_RUN_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$MAASLIN_RUN_ROOT/source"
mode="${1:?prepare, run or collect}"
if [[ -f "$MAASLIN_RUN_ROOT/identity.json" ]]; then
    python3 - "$MAASLIN_RUN_ROOT/identity.json" "$ANALYSIS_SIF" <<'PY'
import json,sys,os
from pathlib import Path
state=json.loads(Path(sys.argv[1]).read_text())
if os.path.abspath(sys.argv[2])!=state['image_path']:
    raise SystemExit('Use the original recorded analysis image for this run')
PY
fi
runtime=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PWD"
    --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,REQUIRE_MAASLIN_BATCH_BACKEND=1
    "$ANALYSIS_SIF")
if [[ "$mode" == prepare ]]; then
    # If the operational gate fails, publish a clear status; dependent fits stay held.
    trap 'mkdir -p "$MAASLIN_RUN_ROOT/PREP_FAILURE"; printf "Preparation failed; see da3_maaslin_prepare.%s.err\n" "$SLURM_JOB_ID" > "$MAASLIN_RUN_ROOT/PREP_FAILURE/README.txt"' ERR
    if [[ "${MAASLIN_RESUME:-0}" == 0 ]]; then
        : "${DA3_BATCH_PLAN:?}"
        sha256sum "$ANALYSIS_SIF" > "$MAASLIN_RUN_ROOT/image.sha256"
        "${runtime[@]}" python3 -B analysis_v2/scripts/maaslin_production.py prepare \
            --source "$DA3_BATCH_PLAN" --root "$MAASLIN_RUN_ROOT" --repo "$PWD" --image "$ANALYSIS_SIF"
    fi
    sha256sum -c --quiet "$MAASLIN_RUN_ROOT/image.sha256"
    (cd "$MAASLIN_RUN_ROOT" && sha256sum -c --quiet identity.sha256)
    "${runtime[@]}" Rscript analysis_v2/tests/test_maaslin_production.R
    echo '[PASS] Operational gate: actual backend and compact output equivalence'
elif [[ "$mode" == run ]]; then
    sha256sum -c --quiet "$MAASLIN_RUN_ROOT/image.sha256"
    index=$(( ${MAASLIN_BATCH_OFFSET:?} + ${SLURM_ARRAY_TASK_ID:?} ))
    "${runtime[@]}" python3 -B analysis_v2/scripts/maaslin_production.py run \
        --root "$MAASLIN_RUN_ROOT" --repo "$PWD" --index "$index"
elif [[ "$mode" == collect ]]; then
    sha256sum -c --quiet "$MAASLIN_RUN_ROOT/image.sha256"
    "${runtime[@]}" python3 -B analysis_v2/scripts/maaslin_production.py collect \
        --root "$MAASLIN_RUN_ROOT" --repo "$PWD" --out "${MAASLIN_REPORT_ROOT:?}"
else
    echo "Unknown mode: $mode" >&2; exit 1
fi
