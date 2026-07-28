#!/usr/bin/env bash
# Purpose: submit Kraken2/Bracken and MetaPhlAn 4 jobs per sample
#          in batches of N samples, wait for each batch to finish, then postprocess.

set -euo pipefail
IFS=$'\n\t'

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${ROOT:-$(cd "$HERE/.." && pwd)}"
ROOT="$(cd "$ROOT" && pwd -P)"

# --- Load configs ---
GLOBAL_ENV="${GLOBAL_ENV:-$ROOT/config/global.env}"
[[ -s "$GLOBAL_ENV" ]] || {
  echo "[error] Missing global configuration: $GLOBAL_ENV" >&2
  echo "        Copy config/global.env.example to config/global.env and edit it." >&2
  exit 1
}
# shellcheck source=/dev/null
source "$GLOBAL_ENV"

: "${DATASET_ENV:?Set DATASET_ENV to a dataset environment file}"
[[ -s "$DATASET_ENV" ]] || {
  echo "[error] Missing dataset configuration: $DATASET_ENV" >&2
  exit 1
}
echo "[info] Using dataset env: $DATASET_ENV"
# shellcheck source=/dev/null
source "$DATASET_ENV"

: "${DATASET:?DATASET not set}"
: "${SAMPLES:?SAMPLES not set (bash array)}"
: "${DATA_DIR:?DATA_DIR not set}"
: "${PAIRED:?PAIRED not set (true|false)}"
: "${SIF:?SIF not set}"
: "${BIND:?BIND not set}"

# Kraken2+Bracken
USE_KRAKEN="${USE_KRAKEN:-1}"
if [[ "$USE_KRAKEN" == "1" ]]; then
  : "${K2_DB:?K2_DB not set (Kraken2 DB dir)}"
  READ_LEN="${READ_LEN:-150}"
fi
K2_THREADS="${K2_THREADS:-16}"
K2_MEM="${K2_MEM:-64}"

# MetaPhlAn4
USE_METAPHLAN="${USE_METAPHLAN:-1}"
if [[ "$USE_METAPHLAN" == "1" ]]; then
  : "${MPA_DB:?MPA_DB not set (MetaPhlAn4 marker DB)}"
fi
MPA_THREADS="${MPA_THREADS:-8}"
MPA_MEM="${MPA_MEM:-16}"

SLURM_PARTITION="${SLURM_PARTITION:-}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"
SLURM_QCHECK_SECS="${SLURM_QCHECK_SECS:-20}"

# ---- NEW: batch size ----
BATCH_SIZE="${BATCH_SIZE:-5}"

ts(){ date +"%Y-%m-%dT%H:%M:%S%z"; }

echo "$(ts) [info] Dataset=$DATASET | Samples=${#SAMPLES[@]} | BatchSize=$BATCH_SIZE"
echo "$(ts) [info] Tools: kraken2_bracken=${USE_KRAKEN} metaphlan4=${USE_METAPHLAN}"

# ---------- Job utils ----------
declare -a SUBMITTED=()      # "<jobid>:<sample>:<tool>"
declare -A JOBSTATE=()       # jid -> state
declare -a OK_JOBS=()        # "jid:sample:tool"
declare -a BAD_JOBS=()       # "jid:sample:tool:state"
declare -a ALL_BAD=()        # accumulate across batches

_submit_job() {
  local tool="$1" sbatch_script="$2" sample="$3" threads="$4" mem_gb="$5"
  local env_args=(
    env
    "ROOT=$ROOT" "DATASET=$DATASET" "SAMPLE=$sample"
    "DATA_DIR=$DATA_DIR" "PAIRED=$PAIRED" "SIF=$SIF" "BIND=$BIND"
  )
  [[ -n "${R1:-}" ]] && env_args+=("R1=$R1")
  [[ -n "${R2:-}" ]] && env_args+=("R2=$R2")
  case "$tool" in
    kraken2_bracken)
      env_args+=("K2_DB=$K2_DB" "READ_LEN=$READ_LEN" "K2_THREADS=$threads" "K2_MEM=$mem_gb")
      ;;
    metaphlan4)
      env_args+=("MPA_DB=$MPA_DB" "MPA_THREADS=$threads" "MPA_MEM=$mem_gb")
      ;;
  esac
  local cmd=( sbatch -D "$ROOT" -c "$threads" --mem="${mem_gb}"G )
  [[ -n "$SLURM_PARTITION" ]] && cmd+=( -p "$SLURM_PARTITION" )
  [[ -n "$SLURM_ACCOUNT"   ]] && cmd+=( -A "$SLURM_ACCOUNT" )
  local out; out="$("${env_args[@]}" "${cmd[@]}" --export=ALL "$sbatch_script")"
  local jid; jid="$(awk '{print $NF}' <<<"$out")"
  [[ "$jid" =~ ^[0-9]+$ ]] || { echo "$(ts) [error] Could not parse JobID: $out" >&2; exit 1; }
  echo "$(ts) [submit] $tool | $sample -> JobID $jid"
  SUBMITTED+=("${jid}:${sample}:${tool}")
}

_job_state() {
  local jid="$1"
  # prefer sacct for finished jobs; squeue for running/pending
  local s; s="$(sacct -j "$jid" --format=State%20 -n | head -n1 | awk '{print $1}')" || true
  if [[ -z "$s" ]]; then s="$(squeue -j "$jid" -h -o "%T" || true)"; fi
  echo "$s"
}

_wait_all() {
  echo "$(ts) [info] Waiting for ${#SUBMITTED[@]} job(s) in this batch…"
  local remaining=${#SUBMITTED[@]}

  # init states
  for triple in "${SUBMITTED[@]}"; do
    IFS=':' read -r jid _ _ <<<"$triple"
    JOBSTATE["$jid"]="PENDING"
  done

  while (( remaining > 0 )); do
    sleep "$SLURM_QCHECK_SECS"
    for triple in "${SUBMITTED[@]}"; do
      IFS=':' read -r jid sample tool <<<"$triple"
      case "${JOBSTATE[$jid]}" in
        COMPLETED|COMPLETED+|CD|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|OOM|PREEMPTED|NODE_FAIL|NF) continue ;;
      esac
      local state; state="$(_job_state "$jid")"; [[ -z "$state" ]] && state="UNKNOWN"
      JOBSTATE["$jid"]="$state"
      case "$state" in
        COMPLETED|COMPLETED+|CD)
          OK_JOBS+=("${jid}:${sample}:${tool}")
          ((remaining--))
          ;;
        FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|OOM|PREEMPTED|NODE_FAIL|NF)
          BAD_JOBS+=("${jid}:${sample}:${tool}:${state}")
          ((remaining--))
          ;;
        # else still running
      esac
    done
  done

  echo "$(ts) [info] Batch finished: ${#OK_JOBS[@]} ok, ${#BAD_JOBS[@]} failed."
  if ((${#BAD_JOBS[@]} > 0)); then
    echo "$(ts) [summary] Failures in this batch:"
    for x in "${BAD_JOBS[@]}"; do
      IFS=':' read -r jid sample tool state <<<"$x"
      echo "  - $tool | $sample | Job $jid | state=$state"
    done
  fi
}

# ---------- MAIN: run in batches ----------
total=${#SAMPLES[@]}
for ((i=0; i<total; i+=BATCH_SIZE)); do
  # slice the next batch
  mapfile -t BATCH < <(printf "%s\n" "${SAMPLES[@]:i:BATCH_SIZE}")

  echo "$(ts) [batch] $((i+1))..$((i+${#BATCH[@]})) / $total samples: ${BATCH[*]}"

  # clear per-batch trackers
  SUBMITTED=(); JOBSTATE=(); OK_JOBS=(); BAD_JOBS=()

  # ---- submit only this batch ----
  if [[ "$USE_KRAKEN" == "1" ]]; then
    for SAMPLE in "${BATCH[@]}"; do
      _submit_job "kraken2_bracken" "$ROOT/workflows/kraken2_bracken/classify_bracken.sbatch" "$SAMPLE" "$K2_THREADS" "$K2_MEM"
    done
  fi
  if [[ "$USE_METAPHLAN" == "1" ]]; then
    for SAMPLE in "${BATCH[@]}"; do
      _submit_job "metaphlan4" "$ROOT/workflows/metaphlan4/profile.sbatch" "$SAMPLE" "$MPA_THREADS" "$MPA_MEM"
    done
  fi

  # # ---- wait for the batch (don’t fail-fast) ----
  _wait_all

  # ---- postprocess only successful jobs in this batch ----
  declare -A OK_BY_TOOL=()
  for triple in "${OK_JOBS[@]}"; do
    IFS=':' read -r _jid sample tool <<<"$triple"
    OK_BY_TOOL["$tool"]+="$sample "
  done

  if [[ "$USE_KRAKEN" == "1" ]]; then
    for SAMPLE in ${OK_BY_TOOL["kraken2_bracken"]-}; do
      echo "$(ts) [postprocess] kraken2_bracken | $SAMPLE"
      ROOT="$ROOT" DATASET="$DATASET" SAMPLE="$SAMPLE" \
        bash "$ROOT/workflows/kraken2_bracken/postprocess_local.sh"
    done
  fi
  if [[ "$USE_METAPHLAN" == "1" ]]; then
    for SAMPLE in ${OK_BY_TOOL["metaphlan4"]-}; do
      echo "$(ts) [postprocess] metaphlan4 | $SAMPLE"
      ROOT="$ROOT" DATASET="$DATASET" SAMPLE="$SAMPLE" \
        bash "$ROOT/workflows/metaphlan4/postprocess_local.sh"
    done
  fi

  # accumulate failures across batches (for final nonzero exit)
  if ((${#BAD_JOBS[@]} > 0)); then
    ALL_BAD+=("${BAD_JOBS[@]}")
  fi

  echo "$(ts) [batch] done."
done

# ---- Final summary & exit code ----
if ((${#ALL_BAD[@]} > 0)); then
  echo "$(ts) [final] Some jobs failed overall:"
  for x in "${ALL_BAD[@]}"; do
    IFS=':' read -r jid sample tool state <<<"$x"
    echo "  - $tool | $sample | Job $jid | state=$state"
  done
  echo "$(ts) [done] Outputs in $ROOT/results/$DATASET/{kraken2_bracken,metaphlan4}"
  exit 2
fi

echo "$(ts) [final] All batches completed successfully."
echo "$(ts) [done] Outputs in $ROOT/results/$DATASET/{kraken2_bracken,metaphlan4}"
