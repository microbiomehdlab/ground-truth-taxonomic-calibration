#!/usr/bin/env bash
#SBATCH --job-name=da3_canary
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=da3_canary.%A_%a.out
#SBATCH --error=da3_canary.%A_%a.err
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA3_CANARY_PLAN:?}" "${DA3_CANARY_RESULTS:?}"
cd "$PROJECT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" \
    --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 \
    "$ANALYSIS_SIF" python3 -B analysis_v2/scripts/plan_da3_canary.py assemble \
    --plan "$DA3_CANARY_PLAN" \
    --outdir "$DA3_CANARY_RESULTS/canary_$(printf '%04d' "$SLURM_ARRAY_TASK_ID")" \
    --index "$SLURM_ARRAY_TASK_ID" --fit
