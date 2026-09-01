#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  echo "Usage: $0 --manifest independent.tsv [--max-concurrent N] [--array-indices SPEC] [--delete-verified-inputs]" >&2
}
manifest=""; max_concurrent=1; array_indices=""; delete_inputs=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) manifest="$2"; shift 2 ;;
    --max-concurrent) max_concurrent="$2"; shift 2 ;;
    --array-indices) array_indices="$2"; shift 2 ;;
    --delete-verified-inputs) delete_inputs=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$manifest" ]] || { usage; exit 2; }
[[ "$max_concurrent" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid concurrency" >&2; exit 2; }
PROJECT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}"
: "${YACHIDA_ENV:?export YACHIDA_ENV}"
source "$YACHIDA_ENV"
: "${ASSEMBLY_SENSITIVITY_ROOT:?set in YACHIDA_ENV}"
manifest="$(cd "$(dirname "$manifest")" && pwd -P)/$(basename "$manifest")"
n=$(( $(wc -l < "$manifest") - 1 ))
[[ "$n" -eq 30 ]] || { echo "[ERROR] Frozen assembly-sensitivity manifest must contain 30 samples; observed $n" >&2; exit 1; }
logs="$ASSEMBLY_SENSITIVITY_ROOT/logs"
mkdir -p "$logs"
array="${array_indices:-1-$n}"
[[ "$array" == *%* ]] || array="${array}%${max_concurrent}"
jid="$(env PROJECT="$PROJECT" YACHIDA_ENV="$YACHIDA_ENV" BATCH_MANIFEST="$manifest" DELETE_INPUTS="$delete_inputs" \
  sbatch --parsable --array="$array" --chdir="$PROJECT" \
    --output="$logs/assembly_sensitivity.%A_%a.out" --error="$logs/assembly_sensitivity.%A_%a.err" \
    --export=ALL "$PROJECT/run_yachida_assembly_sensitivity_sample.sbatch")"
echo "[SUBMITTED] experiment=assembly_sensitivity job=$jid array=$array delete_verified_inputs=$delete_inputs"
echo "[INFO] Logs: $logs"
