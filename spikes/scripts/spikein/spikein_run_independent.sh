#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  spikein_run_independent.sh --env spikein.env [--labels Fnuc,Pmic,Pana]

Options:
  --env FILE         spikein environment file
  --labels CSV       optional subset of labels from SPIKE_PANEL
                     example: --labels Fnuc,Pmic,Pana

Notes:
  - submits one array job per selected taxon
  - each array task processes one sample and loops over all FRACTIONS
  - concurrency is controlled by MAX_CONCURRENT (default 5)
EOF
}

ENV=""
LABELS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV="$2"; shift 2 ;;
    --labels) LABELS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown arg: $1" >&2
      usage
      exit 2
      ;;
  esac
done

[[ -n "$ENV" ]] || { usage; exit 2; }
[[ -f "$ENV" ]] || { echo "[ERROR] Missing env file: $ENV" >&2; exit 1; }

# shellcheck disable=SC1090
source "$ENV"

: "${IMG:?set IMG in env}"
: "${WORK:?set WORK in env}"
: "${SAMPLES_TSV:?set SAMPLES_TSV in env}"
: "${SPIKE_PANEL:?set SPIKE_PANEL in env}"

[[ -f "$SAMPLES_TSV" ]] || { echo "[ERROR] Missing SAMPLES_TSV: $SAMPLES_TSV" >&2; exit 1; }
[[ -f "$SPIKE_PANEL" ]] || { echo "[ERROR] Missing SPIKE_PANEL: $SPIKE_PANEL" >&2; exit 1; }

command -v sbatch >/dev/null 2>&1 || { echo "[ERROR] sbatch not found" >&2; exit 1; }

POOLS_DIR="${POOLS_DIR:-$WORK/pools}"
INDEPENDENT_FRACTIONS="${INDEPENDENT_FRACTIONS:-${FRACTIONS:-0.0001,0.0005,0.001,0.005,0.01,0.05}}"
FRACTIONS="$INDEPENDENT_FRACTIONS"
SEED_BASE="${SEED_BASE:-13}"
MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"
MAX_CONCURRENT="${MAX_CONCURRENT:-5}"
BIND="${BIND:-}"

export FRACTIONS SEED_BASE MAMBA_ENV
echo "[INFO] Independent final fractions: $FRACTIONS"

count_reads() {
  local fq="$1"
  if [[ "$fq" == *.gz ]]; then
    gzip -cd -- "$fq" | awk 'END{print int(NR/4)}'
  else
    awk 'END{print int(NR/4)}' "$fq"
  fi
}

N=$(( $(wc -l < "$SAMPLES_TSV") - 1 ))
(( N > 0 )) || { echo "[ERROR] No samples in $SAMPLES_TSV" >&2; exit 1; }

OUTROOT="$WORK/independent"
mkdir -p "$OUTROOT"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEED_HELPER="${SEED_HELPER:-$SCRIPT_DIR/stable_seed.py}"
[[ -s "$SEED_HELPER" ]] || { echo "[ERROR] Missing seed helper: $SEED_HELPER" >&2; exit 1; }
export SEED_HELPER

cp -f "$ENV" "$OUTROOT/spikein.env.used"

declare -A WANT=()
declare -A SEEN=()

if [[ -n "$LABELS" ]]; then
  IFS=',' read -r -a requested <<< "$LABELS"
  for lab in "${requested[@]}"; do
    lab="$(printf '%s' "$lab" | xargs)"
    [[ -n "$lab" ]] || continue
    WANT["$lab"]=1
  done
fi

PANEL_USED="$OUTROOT/panel.used.tsv"
echo -e "label\ttaxon_name\tassembly\tfasta\tweight\turl" > "$PANEL_USED"

submitted=0

while IFS=$'\t' read -r label taxon_name assembly fasta weight url; do
  [[ "$label" == "label" ]] && continue
  [[ -n "${label:-}" ]] || continue
  [[ "$label" =~ ^# ]] && continue

  if (( ${#WANT[@]} > 0 )); then
    [[ -n "${WANT[$label]+x}" ]] || continue
  fi

  SEEN["$label"]=1

  p1="$POOLS_DIR/${label}.pool_1.fq"
  p2="$POOLS_DIR/${label}.pool_2.fq"
  [[ -s "$p1" && -s "$p2" ]] || {
    echo "[ERROR] Missing/empty pool for $label: $p1 / $p2" >&2
    exit 1
  }

  pool1_n="$(count_reads "$p1")"
  pool2_n="$(count_reads "$p2")"
  [[ "$pool1_n" == "$pool2_n" ]] || {
    echo "[ERROR] Pool mate mismatch for $label: POOL1=$pool1_n POOL2=$pool2_n" >&2
    exit 1
  }
  (( pool1_n > 0 )) || {
    echo "[ERROR] Zero-size pool for $label" >&2
    exit 1
  }

  echo -e "${label}\t${taxon_name}\t${assembly}\t${fasta}\t${weight}\t${url}" >> "$PANEL_USED"

  out="$OUTROOT/$label"
  mkdir -p "$out/logs"
  printf '%s\n' "$label" > "$out/label.used.txt"
  cp -f "$ENV" "$out/spikein.env.used"

  jid=$(
    env \
      IMG="$IMG" WORK="$out" SAMPLES_TSV="$SAMPLES_TSV" \
      POOL1="$p1" POOL2="$p2" POOL_SIZE="$pool1_n" LABEL="$label" BIND="$BIND" SEED_HELPER="$SEED_HELPER" \
    sbatch --parsable \
      --array=1-"$N"%${MAX_CONCURRENT} \
      --chdir="$out" \
      --output="$out/logs/${label}.spike.%A_%a.out" \
      --error="$out/logs/${label}.spike.%A_%a.err" \
      --export=ALL \
      "$SCRIPT_DIR/spike_one_taxon_array.sbatch"
  )

  echo "[SUBMITTED] label=$label array_job=$jid max_concurrent=$MAX_CONCURRENT pool_pairs=$pool1_n"
  submitted=$((submitted + 1))
done < "$SPIKE_PANEL"

if (( ${#WANT[@]} > 0 )); then
  missing=0
  for lab in "${!WANT[@]}"; do
    if [[ -z "${SEEN[$lab]+x}" ]]; then
      echo "[ERROR] Requested label not found in SPIKE_PANEL: $lab" >&2
      missing=1
    fi
  done
  [[ "$missing" -eq 0 ]] || exit 1
fi

(( submitted > 0 )) || {
  echo "[ERROR] No taxa selected for submission" >&2
  exit 1
}

echo "[OK] Submitted $submitted independent taxon job(s)"
echo "[INFO] Panel used: $PANEL_USED"
echo "[INFO] Output root: $OUTROOT"
