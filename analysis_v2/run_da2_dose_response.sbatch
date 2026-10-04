#!/bin/bash
#SBATCH --job-name=da2_response
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=da2_response.%A_%a.out
#SBATCH --error=da2_response.%A_%a.err
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA2_ROOT:?}"
cd "$DA2_ROOT/source"
sha256sum -c --quiet SOURCE_SHA256SUMS
sha256sum -c --quiet "$DA2_ROOT/image.sha256"
args=("$1" --root "$DA2_ROOT" --repo "$PWD")
if [[ "$1" == prepare ]]; then
    args+=(--inventory "${DA_INVENTORY_ROOT:?}" --bracken "${BRACKEN_RECOVERY_ROOT:?}" --metaphlan "${METAPHLAN_RECOVERY_ROOT:?}")
elif [[ "$1" == run ]]; then
    args+=(--index "${SLURM_ARRAY_TASK_ID:?}")
fi
apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1 \
    --bind /mnt/nfs:/mnt/nfs --pwd "$PWD" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/da2_dose_response.py "${args[@]}"
sha256sum -c --quiet SOURCE_SHA256SUMS
