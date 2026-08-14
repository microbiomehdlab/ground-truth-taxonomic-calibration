#!/usr/bin/env bash
# One-time fail-closed gate before submitting Yachida production batches.
set -euo pipefail
IFS=$'\n\t'

ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
YACHIDA_ENV="${YACHIDA_ENV:-$ROOT/config/yachida.env}"
[[ -s "$YACHIDA_ENV" ]] || { echo "[ERROR] Missing site configuration: $YACHIDA_ENV" >&2; exit 1; }
source "$YACHIDA_ENV"

: "${YACHIDA_POOLS_DIR:?set YACHIDA_POOLS_DIR}"
: "${BRACKEN_THRESHOLD:?set BRACKEN_THRESHOLD in YACHIDA_ENV}"
[[ "$BRACKEN_THRESHOLD" == 10 ]] || {
  echo "[ERROR] The frozen Yachida protocol requires BRACKEN_THRESHOLD=10; observed $BRACKEN_THRESHOLD" >&2
  exit 1
}
SPIKE_ENV="${SPIKE_ENV:-$ROOT/work/yachida_67x3/spikein.env}"
[[ -s "$SPIKE_ENV" ]] || { echo "[ERROR] Missing spike environment: $SPIKE_ENV" >&2; exit 1; }
source "$SPIKE_ENV"

PILOT_MANIFEST="${PILOT_MANIFEST:-$ROOT/work/yachida_67x3/metadata/pilot_batched.tsv}"
INDEPENDENT_MANIFEST="${INDEPENDENT_MANIFEST:-$ROOT/work/yachida_67x3/metadata/independent_10_per_condition.tsv}"
LOCAL_PANEL="${LOCAL_PANEL:-$ROOT/spikes/spike_panel.tsv}"
# The production runner derives its equal-weight community from this panel.
COMMUNITY_TSV="${COMMUNITY_TSV:-$LOCAL_PANEL}"
AUDIT_DIR="${PREPRODUCTION_AUDIT_DIR:-$ROOT/work/yachida_67x3/preproduction_audit}"
VERIFY_POOL_CONTENTS="${VERIFY_POOL_CONTENTS:-1}"

for path in "$PILOT_MANIFEST" "$INDEPENDENT_MANIFEST" "$LOCAL_PANEL" "$COMMUNITY_TSV"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing audit input: $path" >&2; exit 1; }
done
arguments=(
  --pilot-manifest "$PILOT_MANIFEST"
  --independent-manifest "$INDEPENDENT_MANIFEST"
  --panel "$LOCAL_PANEL"
  --community "$COMMUNITY_TSV"
  --pools-dir "$YACHIDA_POOLS_DIR"
  --allocator "$ROOT/spikes/scripts/spikein/allocate_community_reads.py"
  --paired-fastq-validator "$ROOT/scripts/validate_paired_fastq.py"
  --output-dir "$AUDIT_DIR"
)
[[ "$VERIFY_POOL_CONTENTS" == 1 ]] && arguments+=(--verify-pool-contents)
python3 "$ROOT/datasets/yachida/preproduction_audit.py" "${arguments[@]}"
