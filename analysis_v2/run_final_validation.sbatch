#!/bin/bash
#SBATCH --job-name=da_validation
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=da_validation.%A_%a.out
#SBATCH --error=da_validation.%A_%a.err
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA2_PAIRED_RESULTS:?}" "${DA_INVENTORY_ROOT:?}" "${FINAL_VALIDATION_ROOT:?}" "${SLURM_ARRAY_TASK_ID:?}"
cd "$PROJECT"
sha256sum -c --quiet "$FINAL_VALIDATION_ROOT/submission.sha256"
apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1 \
  --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
  python3 -B analysis_v2/scripts/run_final_validation_task.py --repo "$PROJECT" \
  --output "$FINAL_VALIDATION_ROOT" --paired "$DA2_PAIRED_RESULTS" \
  --inventory "$DA_INVENTORY_ROOT" --index "$SLURM_ARRAY_TASK_ID"
sha256sum -c --quiet "$FINAL_VALIDATION_ROOT/submission.sha256"
