#!/usr/bin/env bash
#SBATCH --job-name=yachida_final_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=yachida_final_audit.%j.out
#SBATCH --error=yachida_final_audit.%j.err

set -euo pipefail
ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/datasets/yachida/audit_production.py" ]] || {
  echo "[ERROR] Submit from the project root or export PROJECT" >&2
  exit 1
}
cd "$ROOT"
: "${YACHIDA_ENV:?export YACHIDA_ENV}"
source "$YACHIDA_ENV"
: "${YACHIDA_SCRATCH_ROOT:?set in YACHIDA_ENV}"
: "${YACHIDA_STATE_DIR:?set in YACHIDA_ENV}"

python3 "$ROOT/datasets/yachida/audit_production.py" \
  --manifest "$ROOT/work/yachida_67x3/metadata/pilot_batched.tsv" \
  --independent-manifest "$ROOT/work/yachida_67x3/metadata/independent_10_per_condition.tsv" \
  --scratch-root "$YACHIDA_SCRATCH_ROOT" \
  --results-root "$PERSISTENT_RESULTS_ROOT" \
  --state-dir "$YACHIDA_STATE_DIR" \
  --expected-samples 201
