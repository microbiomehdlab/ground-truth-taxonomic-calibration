#!/usr/bin/env bash
#SBATCH --job-name=reference_scale
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=reference_scale.%A_%a.out
#SBATCH --error=reference_scale.%A_%a.err
set -euo pipefail
: "${REFERENCE_SCALE_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$REFERENCE_SCALE_ROOT/source"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --bind /mnt/beegfs:/mnt/beegfs --pwd "$PWD" \
    --env TZ=UTC,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 \
    "$ANALYSIS_SIF" python3 -B analysis_v2/scripts/reference_scaling_sensitivity.py run \
    --root "$REFERENCE_SCALE_ROOT" --index "${SLURM_ARRAY_TASK_ID:?}" --image "$ANALYSIS_SIF"
