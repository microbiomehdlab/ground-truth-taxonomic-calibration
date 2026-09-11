#!/usr/bin/env bash
# Submit a rolling CRC pipeline with bounded downloads and 42 sequential lanes.
set -euo pipefail
IFS=$'\n\t'

usage() {
  echo "Usage: $0 --manifest production.tsv [--lanes 42] [--download-concurrent 8] [--exclude NODES] [--delete-verified-inputs]" >&2
}
manifest=""; lanes=42; downloads=8; exclude=""; delete_inputs=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) manifest="$2"; shift 2 ;;
    --lanes) lanes="$2"; shift 2 ;;
    --download-concurrent) downloads="$2"; shift 2 ;;
    --exclude) exclude="$2"; shift 2 ;;
    --delete-verified-inputs) delete_inputs=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -s "$manifest" ]] || { usage; exit 2; }
[[ "$lanes" =~ ^[1-9][0-9]*$ && "$downloads" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] lanes and download concurrency must be positive integers" >&2; exit 2;
}
PROJECT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)}"
: "${CRC_ENV:?export CRC_ENV}"
manifest="$(cd "$(dirname "$manifest")" && pwd -P)/$(basename "$manifest")"
n=$(( $(wc -l < "$manifest") - 1 ))
(( n > 0 )) || { echo "[ERROR] Empty production manifest" >&2; exit 1; }
cohort="$(basename "$manifest" .tsv)"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
logs="${CRC_ROLLING_LOG_ROOT:-$PROJECT/work/crc_production/logs/${cohort}_rolling_$stamp}"
mkdir -p "$logs"
exclude_args=()
[[ -n "$exclude" ]] && exclude_args+=(--exclude="$exclude")

previous_download=""; previous_compute=""; generation=0
printf 'generation\toffset\ttasks\tdownload_job\tcompute_job\n' > "$logs/jobs.tsv"
for (( offset=0; offset<n; offset+=lanes )); do
  generation=$((generation + 1))
  remaining=$((n - offset)); tasks=$lanes
  (( remaining < tasks )) && tasks=$remaining
  dependency_args=()
  if [[ -n "$previous_download" ]]; then
    # Do not overlap download generations. afterany prevents one failed download
    # from globally blocking later lanes; aftercorr still blocks only the matching
    # lane until its previous compute task succeeds.
    dependency_args+=(--dependency="afterany:${previous_download},aftercorr:${previous_compute}")
  fi
  download_job="$(env PROJECT="$PROJECT" CRC_ENV="$CRC_ENV" BATCH_MANIFEST="$manifest" \
    MANIFEST_INDEX_OFFSET="$offset" sbatch --parsable --array="1-${tasks}%${downloads}" \
    --chdir="$PROJECT" "${exclude_args[@]}" "${dependency_args[@]}" \
    --output="$logs/download.g${generation}.%A_%a.out" \
    --error="$logs/download.g${generation}.%A_%a.err" --export=ALL \
    "$PROJECT/run_crc_download_sample.sbatch")"
  compute_job="$(env PROJECT="$PROJECT" CRC_ENV="$CRC_ENV" BATCH_MANIFEST="$manifest" \
    MANIFEST_INDEX_OFFSET="$offset" CRC_STAGED_ONLY=1 DELETE_INPUTS="$delete_inputs" \
    sbatch --parsable --array="1-${tasks}%${lanes}" --dependency="aftercorr:${download_job}" \
    --chdir="$PROJECT" "${exclude_args[@]}" \
    --output="$logs/compute.g${generation}.%A_%a.out" \
    --error="$logs/compute.g${generation}.%A_%a.err" --export=ALL \
    "$PROJECT/run_crc_batch_sample.sbatch")"
  printf '%s\t%s\t%s\t%s\t%s\n' "$generation" "$offset" "$tasks" \
    "$download_job" "$compute_job" >> "$logs/jobs.tsv"
  echo "[SUBMITTED] generation=$generation indices=$((offset + 1))-$((offset + tasks)) download=$download_job compute=$compute_job"
  previous_download="$download_job"; previous_compute="$compute_job"
done
echo "[INFO] Rolling job ledger: $logs/jobs.tsv"
echo "[INFO] At most $downloads downloads and $lanes compute lanes are active"
