#!/usr/bin/env bash
#SBATCH --job-name=reference_scale_prepare
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=reference_scale_prepare.%j.out
#SBATCH --error=reference_scale_prepare.%j.err
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${REFERENCE_SCALE_ROOT:?}" "${ANALYSIS_SIF:?}"
cd "$PROJECT"
python3 -B analysis_v2/scripts/reference_scaling_sensitivity.py prepare \
    --base "$REFERENCE_ROOT" --out "$REFERENCE_SCALE_ROOT" --image "$ANALYSIS_SIF"
