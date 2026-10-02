#!/usr/bin/env bash
#SBATCH --job-name=preprod_check
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=preproduction.%A_%a.out
#SBATCH --error=preproduction.%A_%a.err
set -euo pipefail
: "${CHECK_ROOT:?}" "${ANALYSIS_SIF:?}" "${CHECK_MODE:?}"
repo="$CHECK_ROOT/source"
cd "$repo"
args=(exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$repo" --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 "$ANALYSIS_SIF" python3 -B)
if [[ "$CHECK_MODE" == validation ]]; then
    apptainer "${args[@]}" analysis_v2/scripts/remaining_validation.py run \
      --plan "$CHECK_ROOT/validation_plan" --results "$CHECK_ROOT/validation_results" \
      --repo "$repo" --index "$SLURM_ARRAY_TASK_ID"
elif [[ "$CHECK_MODE" == batch_canary ]]; then
    index=$(( (SLURM_ARRAY_TASK_ID / 8) * 800 + SLURM_ARRAY_TASK_ID % 8 ))
    apptainer "${args[@]}" analysis_v2/scripts/plan_da3_batches.py run \
      --plan "$CHECK_ROOT/draft_da3_plan" --results "$CHECK_ROOT/batch_canary_results" \
      --repo "$repo" --index "$index"
else
    echo "Unknown check mode" >&2; exit 1
fi
