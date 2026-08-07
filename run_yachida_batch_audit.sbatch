#!/usr/bin/env bash
#SBATCH --job-name=yachida_batch_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=yachida_batch_audit.%j.out
#SBATCH --error=yachida_batch_audit.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export YACHIDA_ENV}"
: "${BATCH_MANIFEST:?export BATCH_MANIFEST}"
source "$YACHIDA_ENV"
: "${YACHIDA_SCRATCH_ROOT:?set in YACHIDA_ENV}"
: "${YACHIDA_STATE_DIR:?set in YACHIDA_ENV}"

python3 "$ROOT/datasets/yachida/audit_batch.py" \
  --manifest "$BATCH_MANIFEST" \
  --scratch-root "$YACHIDA_SCRATCH_ROOT" \
  --state-dir "$YACHIDA_STATE_DIR"
