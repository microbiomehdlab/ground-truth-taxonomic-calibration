#!/usr/bin/env bash
#SBATCH --job-name=yachida_full_sample
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export YACHIDA_ENV}"
: "${BATCH_MANIFEST:?provided by submit_yachida_batch.sh}"
source "$YACHIDA_ENV"
: "${YACHIDA_SCRATCH_ROOT:?set in YACHIDA_ENV}"
: "${YACHIDA_STATE_DIR:?set in YACHIDA_ENV}"

line="$(awk -F'\t' -v task="${SLURM_ARRAY_TASK_ID:?}" 'NR > 1 && ++seen == task {print; exit}' "$BATCH_MANIFEST")"
[[ -n "$line" ]] || { echo "[ERROR] No batch row for array task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
sample_id="$(printf '%s\n' "$line" | cut -f1)"
echo "[INFO] Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Batch manifest: $BATCH_MANIFEST"
echo "[INFO] Sample: $sample_id"
echo "[INFO] Delete verified scratch: ${DELETE_INPUTS:-0}"

args=(
  --manifest "$BATCH_MANIFEST"
  --sample-id "$sample_id"
  --scratch-root "$YACHIDA_SCRATCH_ROOT"
  --state-dir "$YACHIDA_STATE_DIR/samples"
  --runner "$ROOT/datasets/yachida/run_full_streaming_sample.sh"
)
[[ "${DELETE_INPUTS:-0}" == 1 ]] && args+=(--delete-inputs-after-verification)
python3 "$ROOT/datasets/yachida/stream_sample.py" "${args[@]}"
