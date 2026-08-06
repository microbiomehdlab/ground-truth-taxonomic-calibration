#!/usr/bin/env bash
#SBATCH --job-name=yachida_spike_smoke
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=yachida_spike_smoke.%j.out
#SBATCH --error=yachida_spike_smoke.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export the absolute path to config/yachida.env}"
: "${SAMPLE_ID:?export the retained preprocessing-smoke sample ID}"

export PROJECT="$ROOT"
export SPIKE_ENV="${SPIKE_ENV:-$ROOT/work/yachida_67x3/spikein.env}"
export SAMPLE_WORK="${SAMPLE_WORK:-$ROOT/work/yachida_67x3/smoke_scratch/$SAMPLE_ID}"
export YACHIDA_SPIKE_SMOKE_ROOT="${YACHIDA_SPIKE_SMOKE_ROOT:-$ROOT/work/yachida_67x3/spike_smoke/$SAMPLE_ID}"

echo "[INFO] Job: ${SLURM_JOB_ID:-manual}"
echo "[INFO] Sample: $SAMPLE_ID"
echo "[INFO] Spike environment: $SPIKE_ENV"
bash "$ROOT/datasets/yachida/run_spike_smoke_test.sh"
