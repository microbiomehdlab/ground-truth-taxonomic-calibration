#!/usr/bin/env bash
#SBATCH --job-name=crc_download
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=1-00:00:00

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${CRC_ENV:?export CRC_ENV}"
: "${BATCH_MANIFEST:?provided by rolling submitter}"
source "$CRC_ENV"

manifest_index=$(( ${MANIFEST_INDEX_OFFSET:-0} + ${SLURM_ARRAY_TASK_ID:?} ))
line="$(awk -F'\t' -v task="$manifest_index" 'NR > 1 && ++seen == task {print; exit}' "$BATCH_MANIFEST")"
[[ -n "$line" ]] || { echo "[ERROR] No manifest row for index $manifest_index" >&2; exit 1; }
sample_id="$(printf '%s\n' "$line" | cut -f1)"
echo "[INFO] Download job: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Manifest index: $manifest_index"
echo "[INFO] Sample: $sample_id"

python3 "$ROOT/datasets/stream_multirun_sample.py" \
  --manifest "$BATCH_MANIFEST" \
  --sample-id "$sample_id" \
  --scratch-root "$CRC_SCRATCH_ROOT" \
  --state-dir "$CRC_STATE_DIR/samples" \
  --runner "$ROOT/datasets/yachida/run_full_streaming_sample.sh" \
  --stage-only
