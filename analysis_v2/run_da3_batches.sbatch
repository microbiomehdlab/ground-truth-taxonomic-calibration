#!/usr/bin/env bash
#SBATCH --job-name=da3_batches
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=da3_batches.%A_%a.out
#SBATCH --error=da3_batches.%A_%a.err
set -euo pipefail
: "${DA3_BATCH_SOURCE:?}" "${DA3_BATCH_PLAN:?}" "${DA3_BATCH_RESULTS:?}" "${ANALYSIS_SIF:?}"
index=$(( ${DA3_BATCH_OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
cd "$DA3_BATCH_SOURCE"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$DA3_BATCH_SOURCE" \
 --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 "$ANALYSIS_SIF" \
 python3 -B analysis_v2/scripts/plan_da3_batches.py run --plan "$DA3_BATCH_PLAN" \
 --results "$DA3_BATCH_RESULTS" --repo "$DA3_BATCH_SOURCE" --index "$index"
