#!/usr/bin/env bash
#SBATCH --job-name=clinical_wild
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=clinical_wild.%A_%a.out
#SBATCH --error=clinical_wild.%A_%a.err
set -euo pipefail
: "${WILD_ROOT:?}" "${ANALYSIS_SIF:?}"
sha256sum -c --quiet "$WILD_ROOT/image.sha256"
repo="$WILD_ROOT/source"
cd "$repo"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$repo" \
 --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1 "$ANALYSIS_SIF" \
 python3 -B analysis_v2/scripts/clinical_wild_checks.py run \
 --plan "$WILD_ROOT/plan" --results "$WILD_ROOT/results" --index "$SLURM_ARRAY_TASK_ID" --repo "$repo"
