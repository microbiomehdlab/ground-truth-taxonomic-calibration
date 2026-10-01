#!/bin/bash
#SBATCH --job-name=da_pilot_prepare
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=da_pilot_prepare.%j.out
#SBATCH --error=da_pilot_prepare.%j.err
set -euo pipefail
: "${PROJECT:?}"
: "${ANALYSIS_SIF:?}"
: "${DA_INVENTORY_ROOT:?}"
: "${DA1_CHECKPOINT_ROOT:?}"
: "${DA_PILOT_PLAN:?}"
cd "$PROJECT"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
    python3 -B analysis_v2/scripts/prepare_da_pilot.py --inventory-root "$DA_INVENTORY_ROOT" \
    --da1-root "$DA1_CHECKPOINT_ROOT" --outdir "$DA_PILOT_PLAN"
git rev-parse HEAD > "$DA_PILOT_PLAN/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$DA_PILOT_PLAN/analysis_image.sha256"
(cd "$DA_PILOT_PLAN"; find . -maxdepth 1 -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Pilot plan and image provenance verified: $DA_PILOT_PLAN"
