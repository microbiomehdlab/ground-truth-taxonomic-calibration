#!/usr/bin/env bash
#SBATCH --job-name=build_gtc_images
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=build_gtc_images.%j.out
#SBATCH --error=build_gtc_images.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${UPSTREAM_SIF:?export final upstream SIF path when submitting}"
: "${ANALYSIS_SIF:?export final analysis SIF path when submitting}"

export BUILD_TMPDIR="${BUILD_TMPDIR:-${SLURM_TMPDIR:-/tmp}}"
echo "[INFO] Job: ${SLURM_JOB_ID:-manual}"
echo "[INFO] Host: $(hostname)"
echo "[INFO] Project: $ROOT"
echo "[INFO] Build temporary directory: $BUILD_TMPDIR"
echo "[INFO] Upstream image: $UPSTREAM_SIF"
echo "[INFO] Analysis image: $ANALYSIS_SIF"

bash "$ROOT/build_reproducibility_containers.sh"
