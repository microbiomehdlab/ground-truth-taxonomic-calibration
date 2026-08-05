#!/usr/bin/env bash
#SBATCH --job-name=yachida_base_smoke
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=yachida_base_smoke.%j.out
#SBATCH --error=yachida_base_smoke.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export the absolute path to config/yachida.env}"
: "${SAMPLE_ID:?export the retained smoke-test sample ID}"

export SAMPLE_WORK="${SAMPLE_WORK:-$ROOT/work/yachida_67x3/smoke_scratch/$SAMPLE_ID}"
export YACHIDA_BASELINE_SMOKE_ROOT="${YACHIDA_BASELINE_SMOKE_ROOT:-$ROOT/work/yachida_67x3/baseline_smoke}"

echo "[INFO] Job: ${SLURM_JOB_ID:-manual}"
echo "[INFO] Sample: $SAMPLE_ID"
echo "[INFO] Retained preprocessing directory: $SAMPLE_WORK"
bash "$ROOT/datasets/yachida/run_baseline_profiling_smoke.sh"
