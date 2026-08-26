#!/usr/bin/env bash
# Submit a frozen CRC cohort production manifest as a storage-limited array.
set -euo pipefail
IFS=$'\n\t'

usage() {
  echo "Usage: $0 --manifest production.tsv [--max-concurrent N] [--array-indices SPEC] [--delete-verified-inputs]" >&2
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
PROJECT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)}"
: "${CRC_ENV:?export CRC_ENV}"
manifest="$(cd "$(dirname "$manifest")" && pwd -P)/$(basename "$manifest")"
n=$(( $(wc -l < "$manifest") - 1 ))
(( n > 0 )) || { echo "[ERROR] Empty production manifest" >&2; exit 1; }
cohort="$(basename "$manifest" .tsv)"
logs="$PROJECT/work/crc_production/logs/$cohort"
mkdir -p "$logs"
array="${array_indices:-1-$n}"
[[ "$array" == *%* ]] || array="${array}%${max_concurrent}"
jid="$(env PROJECT="$PROJECT" CRC_ENV="$CRC_ENV" BATCH_MANIFEST="$manifest" DELETE_INPUTS="$delete_inputs" \
  sbatch --parsable --array="$array" --chdir="$PROJECT" \
    --output="$logs/${cohort}.%A_%a.out" --error="$logs/${cohort}.%A_%a.err" \
    --export=ALL "$PROJECT/run_crc_batch_sample.sbatch")"
echo "[SUBMITTED] cohort=$cohort job=$jid array=$array delete_verified_inputs=$delete_inputs"
echo "[INFO] Logs: $logs"
