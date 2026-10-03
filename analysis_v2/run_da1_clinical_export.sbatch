#!/bin/bash
#SBATCH --job-name=da1_clinical
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=da1_clinical.%j.out
#SBATCH --error=da1_clinical.%j.err
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA1_EXPORT_ROOT:?}" "${DA_PILOT_PLAN:?}" "${DA_PILOT_RESULTS:?}"
source_root="$DA1_EXPORT_ROOT/source"
cd "$PROJECT"
sha256sum -c --quiet "$DA1_EXPORT_ROOT/analysis_image.sha256"
(cd "$source_root" && sha256sum -c --quiet SOURCE_SHA256SUMS)
# The original plan also binds the container used for the existing fits.
sha256sum -c --quiet "$DA_PILOT_PLAN/analysis_image.sha256"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$source_root" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/export_da1_clinical.py \
    --plan "$DA_PILOT_PLAN" --results "$DA_PILOT_RESULTS" \
    --outdir "$DA1_EXPORT_ROOT/REPORT" --repo "$source_root"
(cd "$source_root" && sha256sum -c --quiet SOURCE_SHA256SUMS)
sha256sum -c --quiet "$DA1_EXPORT_ROOT/analysis_image.sha256"
(cd "$DA1_EXPORT_ROOT/REPORT" && sha256sum -c --quiet SHA256SUMS)
cat "$DA1_EXPORT_ROOT/REPORT/status.json"
