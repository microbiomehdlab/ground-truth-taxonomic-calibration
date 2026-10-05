#!/usr/bin/env bash
#SBATCH --job-name=reference_prepare
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=reference_prepare.%j.out
#SBATCH --error=reference_prepare.%j.err
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}"
: "${DA_INVENTORY_ROOT:?}" "${GEFF_ROOT:?}" "${DA3_OBSERVED_REPORT:?}"
cd "$PROJECT"
runtime=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --bind /mnt/beegfs:/mnt/beegfs
    --pwd "$PROJECT" --env TZ=UTC,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 "$ANALYSIS_SIF")
"${runtime[@]}" python3 -B analysis_v2/tests/test_reference_response_pilot.py
"${runtime[@]}" Rscript analysis_v2/tests/test_reference_response.R "$PROJECT" --actual
"${runtime[@]}" python3 -B analysis_v2/scripts/reference_response_pilot.py prepare \
    --plan "$DA3_BATCH_PLAN" --inventory "$DA_INVENTORY_ROOT" --geff-root "$GEFF_ROOT" \
    --observed-report "$DA3_OBSERVED_REPORT" --repo "$PROJECT" --image "$ANALYSIS_SIF" \
    --out "$REFERENCE_ROOT" --allocations 20 --arms U --batch-size 10
echo '[PASS] Actual pinned backend and immutable 720-context reference plan'
