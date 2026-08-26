#!/usr/bin/env bash
#SBATCH --job-name=crc_full_sample
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${CRC_ENV:?export CRC_ENV}"
: "${BATCH_MANIFEST:?provided by submit_crc_batch.sh}"
source "$CRC_ENV"
: "${CRC_SCRATCH_ROOT:?set in CRC_ENV}"
: "${CRC_STATE_DIR:?set in CRC_ENV}"
: "${INDEPENDENT_MANIFEST:?set in CRC_ENV}"

line="$(awk -F'\t' -v task="${SLURM_ARRAY_TASK_ID:?}" 'NR > 1 && ++seen == task {print; exit}' "$BATCH_MANIFEST")"
[[ -n "$line" ]] || { echo "[ERROR] No manifest row for array task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
sample_id="$(printf '%s\n' "$line" | cut -f1)"
echo "[INFO] Job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Production manifest: $BATCH_MANIFEST"
echo "[INFO] Sample: $sample_id"
echo "[INFO] Delete verified scratch: ${DELETE_INPUTS:-0}"

# The downstream runner uses this historical variable name; CRC_ENV is the
# authoritative cohort configuration supplied by this entry point.
export YACHIDA_ENV="$CRC_ENV"
args=(
  --manifest "$BATCH_MANIFEST"
  --sample-id "$sample_id"
  --scratch-root "$CRC_SCRATCH_ROOT"
  --state-dir "$CRC_STATE_DIR/samples"
  --runner "$ROOT/datasets/yachida/run_full_streaming_sample.sh"
)
[[ "${DELETE_INPUTS:-0}" == 1 ]] && args+=(--delete-inputs-after-verification)
python3 "$ROOT/datasets/stream_multirun_sample.py" "${args[@]}"
