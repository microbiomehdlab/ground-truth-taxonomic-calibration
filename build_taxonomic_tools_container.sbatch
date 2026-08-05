#!/usr/bin/env bash
#SBATCH --job-name=build_gtc_upstream
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=build_gtc_upstream.%j.out
#SBATCH --error=build_gtc_upstream.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${SIF:?export the final upstream SIF path when submitting}"
export BUILD_TMPDIR="${BUILD_TMPDIR:-${SLURM_TMPDIR:-/tmp}}"

echo "[INFO] Job: ${SLURM_JOB_ID:-manual}"
echo "[INFO] Project: $ROOT"
echo "[INFO] Final image: $SIF"
bash "$ROOT/build_taxonomic_tools_container.sh"
