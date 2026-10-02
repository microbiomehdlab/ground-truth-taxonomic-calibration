#!/bin/bash
#SBATCH --job-name=da_robustness
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=da_robustness.%A_%a.out
#SBATCH --error=da_robustness.%A_%a.err
set -euo pipefail
: "${ROBUSTNESS_SOURCE:?}" "${ROBUSTNESS_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$ROBUSTNESS_SOURCE"
sha256sum -c --quiet "$ROBUSTNESS_ROOT/analysis_image.sha256"
if [[ "${ROBUSTNESS_COLLECT:-0}" == 1 ]];then
  python3 -B analysis_v2/scripts/collect_da_robustness.py "$ROBUSTNESS_ROOT"
else
  : "${FINAL_VALIDATION_ROOT:?}" "${SLURM_ARRAY_TASK_ID:?}"
  apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1 \
    --bind /mnt/nfs:/mnt/nfs --pwd "$ROBUSTNESS_SOURCE" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/run_da_robustness.py --repo "$ROBUSTNESS_SOURCE" \
    --root "$ROBUSTNESS_ROOT" --validation "$FINAL_VALIDATION_ROOT" --index "$SLURM_ARRAY_TASK_ID"
fi
sha256sum -c --quiet "$ROBUSTNESS_ROOT/analysis_image.sha256"
