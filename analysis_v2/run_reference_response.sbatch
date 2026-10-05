#!/usr/bin/env bash
#SBATCH --job-name=reference_response
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=reference_response.%A_%a.out
#SBATCH --error=reference_response.%A_%a.err
set -euo pipefail
: "${REFERENCE_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$REFERENCE_ROOT/source"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --bind /mnt/beegfs:/mnt/beegfs --pwd "$PWD" \
    --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 \
    "$ANALYSIS_SIF" python3 -B analysis_v2/scripts/reference_response_pilot.py \
    run --root "$REFERENCE_ROOT" --index "${SLURM_ARRAY_TASK_ID:?}"
