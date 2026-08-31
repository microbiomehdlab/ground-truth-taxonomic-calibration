#!/usr/bin/env bash
#SBATCH --job-name=yachida_pure_pool
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=yachida_pure_pool.%A_%a.out
#SBATCH --error=yachida_pure_pool.%A_%a.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export YACHIDA_ENV}"
source "$YACHIDA_ENV"
: "${YACHIDA_POOLS_DIR:?set in YACHIDA_ENV}"
: "${PURE_POOL_AUDIT_ROOT:?export PURE_POOL_AUDIT_ROOT}"

task="${SLURM_ARRAY_TASK_ID:?submit as an array}"
POOL_LABEL="$(awk -F'\t' -v row="$((task + 1))" 'NR == row {print $1}' "$YACHIDA_POOLS_DIR/pool_pair_counts.tsv")"
[[ -n "$POOL_LABEL" ]] || { echo "[ERROR] No pool for array task $task" >&2; exit 1; }
export PROJECT="$ROOT" POOL_LABEL
echo "[INFO] Pool: $POOL_LABEL"
echo "[INFO] Audit pairs: ${PURE_POOL_AUDIT_PAIRS:-1000000}"
bash "$ROOT/datasets/yachida/run_pure_pool_audit.sh"
