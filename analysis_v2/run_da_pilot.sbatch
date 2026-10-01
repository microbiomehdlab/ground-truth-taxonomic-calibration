#!/bin/bash
#SBATCH --job-name=da_pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=da_pilot.%A_%a.out
#SBATCH --error=da_pilot.%A_%a.err
set -euo pipefail
: "${PROJECT:?}"
: "${ANALYSIS_SIF:?}"
: "${DA_PILOT_PLAN:?}"
: "${DA_PILOT_RESULTS:?}"
: "${SLURM_ARRAY_TASK_ID:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$PROJECT"
expected_image_sha="$(awk 'NR==1 {print $1}' "$DA_PILOT_PLAN/analysis_image.sha256")"
observed_image_sha="$(sha256sum "$ANALYSIS_SIF" | awk '{print $1}')"
test "$expected_image_sha" = "$observed_image_sha"
apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1 \
    --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/run_da_pilot_task.py --plan "$DA_PILOT_PLAN" \
    --output "$DA_PILOT_RESULTS" --repo "$PROJECT" --index "$SLURM_ARRAY_TASK_ID"
