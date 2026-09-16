#!/usr/bin/env bash
# Automatically admit dependency-blocked CRC download tasks while preserving a
# global download-concurrency limit across rolling generations.
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Usage: supervise_crc_rolling.sh --ledger jobs.tsv [options]

Options:
  --download-concurrent N  Global active/admitted download limit (default: 8)
  --poll-seconds N         Seconds between checks (default: 60)
  --state-dir DIR          Persistent supervisor state directory
  --once                   Perform one admission pass and exit
  --dry-run                Print eligible releases without changing Slurm
EOF
}

ledger=""
limit=8
poll_seconds=60
state_dir=""
once=0
dry_run=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ledger) ledger="$2"; shift 2 ;;
    --download-concurrent) limit="$2"; shift 2 ;;
    --poll-seconds) poll_seconds="$2"; shift 2 ;;
    --state-dir) state_dir="$2"; shift 2 ;;
    --once) once=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -s "$ledger" ]] || { echo "[ERROR] Missing or empty ledger: $ledger" >&2; exit 2; }
[[ "$limit" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid download limit: $limit" >&2; exit 2; }
[[ "$poll_seconds" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid poll interval: $poll_seconds" >&2; exit 2; }

for command in squeue sacct scontrol; do
  command -v "$command" >/dev/null || { echo "[ERROR] Required command not found: $command" >&2; exit 2; }
done

header="$(head -n 1 "$ledger")"
[[ "$header" == $'generation\toffset\ttasks\tdownload_job\tcompute_job' ]] || {
  echo "[ERROR] Unexpected rolling-ledger header: $header" >&2
  exit 2
}

ledger="$(cd "$(dirname "$ledger")" && pwd -P)/$(basename "$ledger")"
state_dir="${state_dir:-$(dirname "$ledger")/supervisor}"
mkdir -p "$state_dir"
release_log="$state_dir/releases.tsv"
event_log="$state_dir/events.log"
lock_file="$state_dir/supervisor.lock"
[[ -e "$release_log" ]] || printf 'timestamp\tdownload_task\tpredecessor_compute\tpredecessor_state\n' > "$release_log"

exec 9>"$lock_file"
if ! flock -n 9; then
  echo "[ERROR] Another supervisor already holds $lock_file" >&2
  exit 1
fi

mapfile -t ledger_rows < <(tail -n +2 "$ledger")
(( ${#ledger_rows[@]} > 0 )) || { echo "[ERROR] Ledger has no job rows" >&2; exit 2; }

declare -a download_jobs compute_jobs
declare -A released_tasks
for row in "${ledger_rows[@]}"; do
  IFS=$'\t' read -r generation offset tasks download_job compute_job <<< "$row"
  [[ "$generation" =~ ^[1-9][0-9]*$ && "$offset" =~ ^[0-9]+$ &&
     "$tasks" =~ ^[1-9][0-9]*$ && "$download_job" =~ ^[0-9]+$ &&
     "$compute_job" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] Malformed ledger row: $row" >&2
    exit 2
  }
  download_jobs+=("$download_job")
  compute_jobs+=("$compute_job")
done

while IFS=$'\t' read -r _ task _; do
  [[ "$task" == "download_task" || -z "$task" ]] && continue
  released_tasks["$task"]=1
done < "$release_log"

all_download_jobs="$(IFS=,; echo "${download_jobs[*]}")"

normalise_state() {
  sed -e 's/+.*//' -e 's/[[:space:]]*$//' <<< "$1"
}

task_state() {
  local task="$1" state
  state="$(sacct -X -n -P -j "$task" --format=State 2>/dev/null | cut -d'|' -f1 | head -n 1 || true)"
  normalise_state "$state"
}

is_nonterminal() {
  case "$1" in
    # An empty accounting result immediately after admission is treated
    # conservatively as active so a restart cannot over-admit downloads.
    ""|PENDING|RUNNING|CONFIGURING|COMPLETING|REQUEUED|RESIZING|SUSPENDED) return 0 ;;
    *) return 1 ;;
  esac
}

running_task_ids() {
  squeue -h -r -j "$all_download_jobs" -t R,CF,CG -o '%i' 2>/dev/null || true
}

count_admitted_releases() {
  local task state count=0
  for task in "${!released_tasks[@]}"; do
    state="$(task_state "$task")"
    is_nonterminal "$state" && count=$((count + 1))
  done
  echo "$count"
}

run_pass() {
  local -A running=()
  local task state download_job previous_compute reason array_job array_task
  local unmanaged_running=0 admitted slots released_now=0

  while IFS= read -r task; do
    [[ -n "$task" ]] && running["$task"]=1
  done < <(running_task_ids)

  for task in "${!running[@]}"; do
    [[ -v "released_tasks[$task]" ]] || unmanaged_running=$((unmanaged_running + 1))
  done
  admitted="$(count_admitted_releases)"
  slots=$((limit - unmanaged_running - admitted))
  (( slots < 0 )) && slots=0

  printf '[%s] running_unmanaged=%d admitted_releases=%d limit=%d slots=%d\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$unmanaged_running" "$admitted" "$limit" "$slots" | tee -a "$event_log"
  (( slots > 0 )) || return 0

  for ((index=1; index<${#download_jobs[@]} && slots>0; index++)); do
    download_job="${download_jobs[index]}"
    previous_compute="${compute_jobs[index-1]}"
    while IFS='|' read -r array_job array_task reason; do
      [[ "$array_task" =~ ^[1-9][0-9]*$ ]] || continue
      [[ "$reason" == Dependency* ]] || continue
      task="${array_job}_${array_task}"
      [[ -v "released_tasks[$task]" ]] && continue
      state="$(task_state "${previous_compute}_${array_task}")"
      [[ "$state" == COMPLETED ]] || continue

      if (( dry_run )); then
        echo "[DRY-RUN] release=$task predecessor=${previous_compute}_${array_task} state=$state"
      else
        scontrol update JobId="$task" Dependency=
        printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
          "$task" "${previous_compute}_${array_task}" "$state" >> "$release_log"
        printf '[RELEASED] %s after verified predecessor %s=%s\n' \
          "$task" "${previous_compute}_${array_task}" "$state" | tee -a "$event_log"
        released_tasks["$task"]=1
      fi
      slots=$((slots - 1))
      released_now=$((released_now + 1))
      (( slots == 0 )) && break
    done < <(squeue -h -r -j "$download_job" -t PD -o '%F|%K|%r' 2>/dev/null || true)
  done
  echo "[INFO] Eligible tasks released this pass: $released_now" | tee -a "$event_log"
}

echo "[INFO] Supervising ledger: $ledger"
echo "[INFO] Persistent release state: $release_log"
echo "[INFO] Global admitted-download limit: $limit"

while true; do
  run_pass
  (( once )) && break
  sleep "$poll_seconds"
done
