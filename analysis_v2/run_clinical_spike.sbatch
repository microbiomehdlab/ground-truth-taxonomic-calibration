#!/usr/bin/env bash
#SBATCH --job-name=clinical_spike
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=clinical_spike.%A_%a.out
#SBATCH --error=clinical_spike.%A_%a.err
set -euo pipefail
: "${CLINICAL_SPIKE_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$CLINICAL_SPIKE_ROOT/source"
sha256sum -c --quiet SOURCE_SHA256SUMS
sha256sum -c --quiet "$CLINICAL_SPIKE_ROOT/image.sha256"
args=("$1" --root "$CLINICAL_SPIKE_ROOT" --repo "$PWD")
if [[ "$1" == prepare ]]; then
    apptainer exec --cleanenv --env REQUIRE_CLINICAL_SPIKE_BACKEND=1 \
        --bind /mnt/nfs:/mnt/nfs --pwd "$PWD" "$ANALYSIS_SIF" \
        Rscript analysis_v2/tests/test_clinical_spike_fit.R
    args+=(--clinical "${DA1_REPORT:?}" --inventory "${DA_INVENTORY_ROOT:?}")
elif [[ "$1" == run ]]; then
    args+=(--index "${SLURM_ARRAY_TASK_ID:?}")
fi
apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1 \
    --bind /mnt/nfs:/mnt/nfs --pwd "$PWD" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/clinical_spike_association.py "${args[@]}"
sha256sum -c --quiet SOURCE_SHA256SUMS
