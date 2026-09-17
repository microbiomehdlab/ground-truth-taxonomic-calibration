#!/usr/bin/env bash
# Admit CRC downloads against a global sample-slot pool instead of fixed lanes.
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Usage: supervise_crc_global_slots.sh --ledger jobs.tsv [options]

Options:
  --inflight-limit N       Maximum admitted samples (default: 42)
  --download-concurrent N  Maximum active/ready downloads (default: 8)
  --external-pipeline D:C  Retry download/compute array pair; repeatable
  --poll-seconds N         Seconds between checks (default: 60)
  --state-dir DIR          Persistent state directory
  --once                   Perform one pass and exit
  --dry-run                Print releases without changing Slurm

An in-flight sample is counted once while downloading, scheduler-ready for
compute, or computing. A completed download whose compute remains dependency-
blocked does not occupy a runnable lane. A completed compute frees its slot.
Only download dependencies are cleared; compute dependencies are never edited.
EOF
}

ledger=""; inflight_limit=42; download_limit=8; poll_seconds=60
state_dir=""; once=0; dry_run=0; external_pipelines=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ledger) ledger="$2"; shift 2 ;;
    --inflight-limit) inflight_limit="$2"; shift 2 ;;
    --download-concurrent) download_limit="$2"; shift 2 ;;
    --external-pipeline) external_pipelines+=("$2"); shift 2 ;;
    --poll-seconds) poll_seconds="$2"; shift 2 ;;
    --state-dir) state_dir="$2"; shift 2 ;;
    --once) once=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -s "$ledger" ]] || { echo "[ERROR] Missing ledger: $ledger" >&2; exit 2; }
[[ "$inflight_limit" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid in-flight limit" >&2; exit 2; }
[[ "$download_limit" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid download limit" >&2; exit 2; }
[[ "$poll_seconds" =~ ^[1-9][0-9]*$ ]] || { echo "[ERROR] Invalid poll interval" >&2; exit 2; }
(( download_limit <= inflight_limit )) || { echo "[ERROR] Download limit exceeds in-flight limit" >&2; exit 2; }
for command in squeue sacct scontrol flock; do
  command -v "$command" >/dev/null || { echo "[ERROR] Missing command: $command" >&2; exit 2; }
done

header="$(head -n 1 "$ledger")"
[[ "$header" == $'generation\toffset\ttasks\tdownload_job\tcompute_job' ]] || {
  echo "[ERROR] Unexpected ledger header: $header" >&2; exit 2;
}
ledger="$(cd "$(dirname "$ledger")" && pwd -P)/$(basename "$ledger")"
state_dir="${state_dir:-$(dirname "$ledger")/global_slots}"
mkdir -p "$state_dir"
release_log="$state_dir/admissions.tsv"
event_log="$state_dir/events.log"
lock_file="$state_dir/supervisor.lock"
[[ -e "$release_log" ]] || printf 'timestamp\tsample_key\tdownload_task\tprevious_dependency\n' > "$release_log"

exec 9>"$lock_file"
flock -n 9 || { echo "[ERROR] Another global-slot supervisor holds $lock_file" >&2; exit 1; }

declare -a download_jobs compute_jobs generations task_counts
declare -A generation_by_download generation_by_compute admitted
while IFS=$'\t' read -r generation offset tasks download_job compute_job; do
  [[ "$generation" == generation ]] && continue
  [[ "$generation" =~ ^[1-9][0-9]*$ && "$tasks" =~ ^[1-9][0-9]*$ &&
     "$download_job" =~ ^[0-9]+$ && "$compute_job" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] Malformed ledger row" >&2; exit 2;
  }
  generations+=("$generation"); task_counts+=("$tasks")
  download_jobs+=("$download_job"); compute_jobs+=("$compute_job")
  generation_by_download["$download_job"]="$generation"
  generation_by_compute["$compute_job"]="$generation"
done < "$ledger"
(( ${#download_jobs[@]} > 0 )) || { echo "[ERROR] Empty ledger" >&2; exit 2; }

while IFS=$'\t' read -r _ key task _; do
  [[ "$key" == sample_key || -z "$key" ]] && continue
  admitted["$key"]="$task"
done < "$release_log"

for pair in "${external_pipelines[@]}"; do
  [[ "$pair" =~ ^[0-9]+:[0-9]+$ ]] || { echo "[ERROR] Invalid external pipeline: $pair" >&2; exit 2; }
done

all_download_jobs="$(IFS=,; echo "${download_jobs[*]}")"
all_compute_jobs="$(IFS=,; echo "${compute_jobs[*]}")"

normalise_state() { sed -e 's/+.*//' -e 's/[[:space:]]*$//' <<< "$1"; }
is_nonterminal() {
  case "$1" in
    PENDING|RUNNING|CONFIGURING|COMPLETING|REQUEUED|RESIZING|SUSPENDED) return 0 ;;
    *) return 1 ;;
  esac
}

# Populate associative maps keyed by displayed array task (JOB_TASK).
load_accounting() {
  local jobs="$1" map_name="$2" line job state
  local -n output="$map_name"
  output=()
  while IFS='|' read -r job state _; do
    [[ "$job" =~ ^[0-9]+_[0-9]+$ ]] || continue
    output["$job"]="$(normalise_state "$state")"
  done < <(sacct -X -n -P -j "$jobs" --format=JobID,State 2>/dev/null || true)
}

load_queue() {
  local jobs="$1" map_name="$2" line array task state reason
  local -n output="$map_name"
  output=()
  while IFS='|' read -r array task state reason; do
    [[ "$array" =~ ^[0-9]+$ && "$task" =~ ^[0-9]+$ ]] || continue
    output["${array}_${task}"]="$state|$reason"
  done < <(squeue -h -r -j "$jobs" -o '%F|%K|%T|%r' 2>/dev/null || true)
}

run_pass() {
  local -A dstate=() cstate=() dqueue=() cqueue=() occupied=() active_download=()
  local -A external_occupied=()
  local i generation tasks djob cjob task key state q reason
  local occupied_count=0 active_download_count=0 free_samples free_downloads admit_now=0

  load_accounting "$all_download_jobs" dstate
  load_accounting "$all_compute_jobs" cstate
  load_queue "$all_download_jobs" dqueue
  load_queue "$all_compute_jobs" cqueue

  for ((i=0; i<${#download_jobs[@]}; i++)); do
    generation="${generations[i]}"; tasks="${task_counts[i]}"
    djob="${download_jobs[i]}"; cjob="${compute_jobs[i]}"
    for ((task=1; task<=tasks; task++)); do
      key="${generation}:${task}"
      state="${cstate[${cjob}_${task}]:-}"
      [[ "$state" == COMPLETED ]] && continue

      q="${cqueue[${cjob}_${task}]:-}"
      reason="${q#*|}"
      if [[ -n "$q" && "$reason" != Dependency* ]]; then occupied["$key"]=1; fi

      state="${dstate[${djob}_${task}]:-}"
      q="${dqueue[${djob}_${task}]:-}"
      reason="${q#*|}"
      if [[ -n "$q" && "$reason" != Dependency* ]]; then
        occupied["$key"]=1
      elif [[ -v "admitted[$key]" ]] && is_nonterminal "$state"; then
        occupied["$key"]=1
      elif (( generation == 1 )) && is_nonterminal "$state"; then
        occupied["$key"]=1
      fi
      if [[ -n "$q" && "$reason" != Dependency* ]]; then active_download["${djob}_${task}"]=1; fi
    done
  done

  # Retry arrays are outside the ledger. Count each array index once across its
  # download->compute transition and include active retry downloads in the cap.
  local pair ed ec ext_jobs ext_task
  for pair in "${external_pipelines[@]}"; do
    ed="${pair%%:*}"; ec="${pair##*:}"; ext_jobs="$ed,$ec"
    local -A estate=() equeue=()
    load_accounting "$ext_jobs" estate
    load_queue "$ext_jobs" equeue
    for ext_task in "${!estate[@]}" "${!equeue[@]}"; do
      [[ "$ext_task" =~ ^(${ed}|${ec})_([0-9]+)$ ]] || continue
      task="${BASH_REMATCH[2]}"; key="external:${ed}:${ec}:${task}"
      [[ "${estate[${ec}_${task}]:-}" == COMPLETED ]] && continue
      q="${equeue[${ec}_${task}]:-}"; reason="${q#*|}"
      # A terminally failed retry download plus DependencyNeverSatisfied
      # compute is not active; it needs a new documented retry.
      if [[ "${estate[${ed}_${task}]:-}" =~ ^(FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED)$ &&
            "$reason" == DependencyNeverSatisfied* ]]; then
        :
      elif is_nonterminal "${estate[${ed}_${task}]:-}" || is_nonterminal "${estate[${ec}_${task}]:-}"; then
        external_occupied["$key"]=1
      fi
      q="${equeue[${ed}_${task}]:-}"; reason="${q#*|}"
      [[ -n "$q" && "$reason" != Dependency* ]] && active_download["${ed}_${task}"]=1
    done
  done

  occupied_count=$(( ${#occupied[@]} + ${#external_occupied[@]} ))
  active_download_count=${#active_download[@]}
  free_samples=$((inflight_limit - occupied_count)); (( free_samples < 0 )) && free_samples=0
  free_downloads=$((download_limit - active_download_count)); (( free_downloads < 0 )) && free_downloads=0
  (( free_samples < free_downloads )) && admit_now=$free_samples || admit_now=$free_downloads

  printf '[%s] occupied=%d external=%d inflight_limit=%d active_downloads=%d download_limit=%d admit=%d\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${#occupied[@]}" "${#external_occupied[@]}" \
    "$inflight_limit" "$active_download_count" "$download_limit" "$admit_now" | tee -a "$event_log"
  (( admit_now > 0 )) || return 0

  # Oldest pending sample first. DependencyNeverSatisfied is admissible here:
  # global tokens replace fixed-lane predecessor gating.
  for ((i=0; i<${#download_jobs[@]} && admit_now>0; i++)); do
    generation="${generations[i]}"; djob="${download_jobs[i]}"
    while IFS='|' read -r array task reason; do
      [[ "$task" =~ ^[0-9]+$ && "$reason" == Dependency* ]] || continue
      key="${generation}:${task}"
      [[ -v "occupied[$key]" || -v "admitted[$key]" ]] && continue
      if (( dry_run )); then
        echo "[DRY-RUN] admit=$key download=${array}_${task} old_reason=$reason"
      else
        scontrol update JobId="${array}_${task}" Dependency=
        printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
          "$key" "${array}_${task}" "$reason" >> "$release_log"
        admitted["$key"]="${array}_${task}"
        echo "[ADMITTED] $key download=${array}_${task} old_reason=$reason" | tee -a "$event_log"
      fi
      admit_now=$((admit_now - 1))
      (( admit_now == 0 )) && break
    done < <(squeue -h -r -j "$djob" -t PD -o '%F|%K|%r' 2>/dev/null || true)
  done
}

echo "[INFO] Global sample slots: $inflight_limit"
echo "[INFO] Global download limit: $download_limit"
echo "[INFO] Ledger: $ledger"
echo "[INFO] State: $state_dir"
while true; do
  run_pass
  (( once )) && break
  sleep "$poll_seconds"
done
