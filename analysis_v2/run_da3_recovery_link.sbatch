#!/usr/bin/env bash
#SBATCH --job-name=da3_recovery_link
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=da3_recovery_link.%j.out
#SBATCH --error=da3_recovery_link.%j.err
set -euo pipefail
: "${RECOVERY_LINK_ROOT:?}" "${RECOVERY_LINK_SOURCE:?}" "${DA3_BATCH_PLAN:?}" "${DA3_RUN_ROOT:?}"
: "${BRACKEN_RECOVERY_ROOT:?}" "${METAPHLAN_RECOVERY_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$RECOVERY_LINK_SOURCE"
sha256sum -c --quiet "$RECOVERY_LINK_ROOT/image.sha256"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$RECOVERY_LINK_SOURCE" \
 "$ANALYSIS_SIF" python3 -B analysis_v2/scripts/link_da3_recovery.py \
 --plan "$DA3_BATCH_PLAN" --results "$DA3_RUN_ROOT/results" \
 --bracken "$BRACKEN_RECOVERY_ROOT" --metaphlan "$METAPHLAN_RECOVERY_ROOT" \
 --out "$RECOVERY_LINK_ROOT/report"
(cd "$RECOVERY_LINK_ROOT/report" && sha256sum -c --quiet SHA256SUMS)
