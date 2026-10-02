#!/bin/bash
#SBATCH --job-name=da3_candidate
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=da3_candidate.%A_%a.out
#SBATCH --error=da3_candidate.%A_%a.err
set -euo pipefail
: "${CANDIDATE_SOURCE:?}" "${CANDIDATE_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$CANDIDATE_SOURCE"
sha256sum -c --quiet "$CANDIDATE_ROOT/analysis_image.sha256"
if [[ "${CANDIDATE_COLLECT:-0}" == 1 ]]; then
  python3 -B analysis_v2/scripts/collect_da3_candidate.py "$CANDIDATE_ROOT"
else
  : "${FINAL_VALIDATION_ROOT:?}" "${DA_PILOT_RESULTS:?}" "${DA_PILOT_PLAN:?}" "${SLURM_ARRAY_TASK_ID:?}"
  apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1 \
    --bind /mnt/nfs:/mnt/nfs --pwd "$CANDIDATE_SOURCE" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/run_da3_candidate.py --repo "$CANDIDATE_SOURCE" \
    --root "$CANDIDATE_ROOT" --validation "$FINAL_VALIDATION_ROOT" \
    --pilot "$DA_PILOT_RESULTS" --plan "$DA_PILOT_PLAN" --index "$SLURM_ARRAY_TASK_ID"
fi
sha256sum -c --quiet "$CANDIDATE_ROOT/analysis_image.sha256"
