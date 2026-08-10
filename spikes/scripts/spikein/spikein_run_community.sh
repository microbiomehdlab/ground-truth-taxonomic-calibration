#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage(){
  cat <<'EOF'
Usage:
  spikein_run_community.sh --env spikein.env --community-label CRCpanel [options]

Options:
  --env FILE              Env file to source
  --community-label STR   Community label (default: CRCpanel)
  --community FILE        Existing community.tsv to use
  --array-indices SPEC    Submit only these sample array indices, e.g. 3,8,12 or 3-12
  --run-fractions LIST    Run only this subset of fractions, e.g. 0.01 or 0.001,0.01
  --time HH:MM:SS         Override sbatch time limit for this submission
  -h, --help              Show help

community.tsv format:
  label<TAB>weight

Notes:
  - Uses a SLURM array over SAMPLES_TSV rows (excluding header)
  - Each array task processes ONE sample and loops over RUN_FRACTIONS
  - FULL_FRACTIONS comes from env FRACTIONS (or default 0.0001,0.001,0.01)
  - RUN_FRACTIONS defaults to FULL_FRACTIONS
  - Seeds derive from stable sample, community, fraction, and taxon identifiers;
    manifest order, array indices, and fraction subsetting cannot change them
  - Concurrency is controlled by MAX_CONCURRENT (default 5): --array=...%MAX_CONCURRENT
EOF
}

ENV=""
COMMUNITY_LABEL="CRCpanel"
COMMUNITY=""
ARRAY_INDICES=""
RUN_FRACTIONS_ARG=""
SBATCH_TIME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV="$2"; shift 2;;
    --community-label) COMMUNITY_LABEL="$2"; shift 2;;
    --community) COMMUNITY="$2"; shift 2;;
    --array-indices) ARRAY_INDICES="$2"; shift 2;;
    --run-fractions) RUN_FRACTIONS_ARG="$2"; shift 2;;
    --time) SBATCH_TIME="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 2 ;;
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

command -v sbatch >/dev/null 2>&1 || { echo "[ERROR] sbatch not found" >&2; exit 1; }

POOLS_DIR="${POOLS_DIR:-$WORK/pools}"
FULL_FRACTIONS="${COMMUNITY_FRACTIONS:-${FRACTIONS:-0.0001,0.0005,0.001,0.005,0.01,0.05,0.10}}"
RUN_FRACTIONS="${RUN_FRACTIONS_ARG:-${RUN_FRACTIONS:-$FULL_FRACTIONS}}"
SEED_BASE="${SEED_BASE:-13}"
MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"
MAX_CONCURRENT="${MAX_CONCURRENT:-5}"
BIND="${BIND:-}"
SAMPLING_MODE="${SAMPLING_MODE:-single_pass}"

# Keep the original full fraction list available for seed reproducibility.
# FRACTIONS is kept for backward compatibility.
export FULL_FRACTIONS RUN_FRACTIONS SEED_BASE MAMBA_ENV
export FRACTIONS="$FULL_FRACTIONS"

N=$(( $(wc -l < "$SAMPLES_TSV") - 1 ))
(( N > 0 )) || { echo "[ERROR] No samples in $SAMPLES_TSV" >&2; exit 1; }

OUTROOT="$WORK/community"
out="$OUTROOT/$COMMUNITY_LABEL"
mkdir -p "$out/logs"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEED_HELPER="${SEED_HELPER:-$SCRIPT_DIR/stable_seed.py}"
[[ -s "$SEED_HELPER" ]] || { echo "[ERROR] Missing seed helper: $SEED_HELPER" >&2; exit 1; }
export SEED_HELPER

# Build community file if not provided
if [[ -z "$COMMUNITY" ]]; then
  COMMUNITY="$out/community.tsv"
  echo -e "label\tweight" > "$COMMUNITY"
  tail -n +2 "$SPIKE_PANEL" | while IFS=$'\t' read -r label taxon_name assembly fasta weight url; do
    [[ -z "${label:-}" ]] && continue
    [[ "$label" =~ ^# ]] && continue
    w="${weight:-1}"
    echo -e "${label}\t${w}" >> "$COMMUNITY"
  done
fi

[[ -f "$COMMUNITY" ]] || { echo "[ERROR] Missing community file: $COMMUNITY" >&2; exit 1; }

# Copy the community used into output dir (provenance)
cp -f "$COMMUNITY" "$out/community.used.tsv"
COMMUNITY="$out/community.used.tsv"

# Validate: weights numeric and pools exist for each label
bad=0
while IFS=$'\t' read -r lab w; do
  [[ -z "${lab:-}" ]] && continue
  [[ "$lab" =~ ^# ]] && continue
  [[ "$lab" == "label" ]] && continue
  w="${w:-1}"

  python3 - <<PY >/dev/null 2>&1
float("$w")
PY
  [[ $? -eq 0 ]] || { echo "[ERROR] Non-numeric weight for $lab: '$w' in $COMMUNITY" >&2; bad=1; }

  p1="$POOLS_DIR/${lab}.pool_1.fq"
  p2="$POOLS_DIR/${lab}.pool_2.fq"
  [[ -s "$p1" && -s "$p2" ]] || { echo "[ERROR] Missing/empty pool for $lab: $p1 / $p2" >&2; bad=1; }
done < "$COMMUNITY"

[[ "$bad" -eq 0 ]] || exit 1

if [[ -n "$ARRAY_INDICES" ]]; then
  if [[ "$ARRAY_INDICES" == *%* ]]; then
    ARRAY_EXPR="$ARRAY_INDICES"
  else
    ARRAY_EXPR="${ARRAY_INDICES}%${MAX_CONCURRENT}"
  fi
else
  ARRAY_EXPR="1-${N}%${MAX_CONCURRENT}"
fi

SBATCH_OPTS=()
[[ -n "$SBATCH_TIME" ]] && SBATCH_OPTS+=("--time=$SBATCH_TIME")

jid=$(
  env \
    IMG="$IMG" WORK="$out" SAMPLES_TSV="$SAMPLES_TSV" \
    COMMUNITY_TSV="$COMMUNITY" POOLS_DIR="$POOLS_DIR" \
    COMMUNITY_LABEL="$COMMUNITY_LABEL" BIND="$BIND" SEED_HELPER="$SEED_HELPER" SAMPLING_MODE="$SAMPLING_MODE" \
  sbatch --parsable \
    "${SBATCH_OPTS[@]}" \
    --array="$ARRAY_EXPR" \
    --chdir="$out" \
    --output="$out/logs/${COMMUNITY_LABEL}.%A_%a.out" \
    --error="$out/logs/${COMMUNITY_LABEL}.%A_%a.err" \
    --export=ALL \
    "$SCRIPT_DIR/spike_community_array.sbatch"
)

echo "[SUBMITTED] community=$COMMUNITY_LABEL array_job=$jid array=$ARRAY_EXPR max_concurrent=$MAX_CONCURRENT"
echo "[INFO] Community file used: $COMMUNITY"
echo "[INFO] FULL_FRACTIONS=$FULL_FRACTIONS"
echo "[INFO] RUN_FRACTIONS=$RUN_FRACTIONS"
echo "[INFO] Sampling mode: $SAMPLING_MODE"
