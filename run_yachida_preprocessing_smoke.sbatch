#!/usr/bin/env bash
#SBATCH --job-name=yachida_prep_smoke
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=yachida_prep_smoke.%j.out
#SBATCH --error=yachida_prep_smoke.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export the absolute path to config/yachida.env}"
: "${SAMPLE_ID:?export one sample ID from the frozen pilot manifest}"

manifest="${YACHIDA_MANIFEST:-$ROOT/work/yachida_67x3/metadata/pilot_batched.tsv}"
scratch_root="${YACHIDA_SMOKE_SCRATCH:-$ROOT/work/yachida_67x3/smoke_scratch}"
state_dir="${YACHIDA_STATE_DIR:-$ROOT/work/yachida_67x3/state/preprocessing_smoke}"

echo "[INFO] Sample: $SAMPLE_ID"
echo "[INFO] Manifest: $manifest"
echo "[INFO] Scratch retained at: $scratch_root"
python3 "$ROOT/datasets/yachida/stream_sample.py" \
  --manifest "$manifest" \
  --sample-id "$SAMPLE_ID" \
  --scratch-root "$scratch_root" \
  --state-dir "$state_dir" \
  --runner "$ROOT/datasets/yachida/run_preprocessing_smoke_test.sh"
