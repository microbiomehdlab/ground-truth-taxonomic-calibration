#!/bin/bash
#SBATCH -J k2b_classify
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o k2b_classify.%x.%j.out
#SBATCH -e k2b_classify.%x.%j.err
#SBATCH --exclude=compute-[1,5,9,13,17,21]

set -euo pipefail
IFS=$'\n\t'
umask 0022
export LC_ALL=C

# -----------------------------
# Resolve repository root
# -----------------------------
if [[ -n "${ROOT:-}" && -d "$ROOT" ]]; then
  :
else
  ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi

# Load globals if present
source "$ROOT/config/global.env" || true

# -----------------------------
# Required environment
# -----------------------------
: "${SAMPLE:?SAMPLE not set}"
: "${DATASET:?DATASET not set}"
: "${DATA_DIR:?DATA_DIR not set}"
: "${SIF:?SIF not set}"
: "${BIND:?BIND not set}"
: "${K2_DB:?K2_DB not set}"

PAIRED="${PAIRED:-false}"

# -----------------------------
# Tunables
# -----------------------------
READ_LEN="${READ_LEN:-150}"
THREADS="${K2_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-10}"
[[ "$BRACKEN_THRESHOLD" =~ ^[1-9][0-9]*$ ]] || {
  echo "[error] BRACKEN_THRESHOLD must be a positive integer" >&2
  exit 1
}
MEM_GB="${K2_MEM:-${SLURM_MEM_PER_NODE:-16}}"

# Prefer direct tool execution in the container.
# Options: auto | direct | micromamba
K2_RUNTIME_MODE="${K2_RUNTIME_MODE:-auto}"
MICROMAMBA_ENV_NAME="${MICROMAMBA_ENV_NAME:-taxonomic_tools}"

# Output retention
# Set to 1 only if downstream code still needs these files.
KEEP_KRAKEN_OUT="${KEEP_KRAKEN_OUT:-0}"
KEEP_KRAKEN_REPORT="${KEEP_KRAKEN_REPORT:-0}"
KEEP_BRACKEN_GENUS="${KEEP_BRACKEN_GENUS:-0}"
KEEP_QC="${KEEP_QC:-0}"

# -----------------------------
# Output paths
# -----------------------------
OUT_ROOT="${OUT_ROOT:-"$ROOT/results"}"
OUT_DIR="$OUT_ROOT/${DATASET}/kraken2_bracken"
mkdir -p "$OUT_DIR"

KRAKEN_REPORT="$OUT_DIR/${SAMPLE}.kraken2.report"
BRACKEN_S="$OUT_DIR/${SAMPLE}.bracken.S.tsv"
BRACKEN_G="$OUT_DIR/${SAMPLE}.bracken.G.tsv"
VERSIONS="$OUT_DIR/${SAMPLE}.versions.txt"
PARAMS_JSON="$OUT_DIR/${SAMPLE}.params.json"
TOTAL_FILE="$OUT_DIR/${SAMPLE}.total_reads.txt"

if [[ "$KEEP_KRAKEN_OUT" == "1" ]]; then
  KRAKEN_OUT="$OUT_DIR/${SAMPLE}.kraken2.out"
else
  KRAKEN_OUT="/dev/null"
fi

# -----------------------------
# Build bind args safely
# -----------------------------
BIND_ARR=()
if [[ -n "${BIND:-}" ]]; then
  if [[ "$BIND" == -* ]]; then
    # e.g. BIND="--bind /mnt/nfs:/mnt/nfs,/scratch:/scratch"
    read -r -a BIND_ARR <<< "$BIND"
  else
    # e.g. BIND="/mnt/nfs:/mnt/nfs,/scratch:/scratch"
    BIND_ARR=(-B "$BIND")
  fi
fi

# Extra environment passed into container
APPT_ENV_ARGS=()

# -----------------------------
# Helpers
# -----------------------------
srun_exec() {
  srun --ntasks=1 --cpus-per-task="$THREADS" --cpu-bind=cores --exact "$@"
}

srun_shell() {
  srun --ntasks=1 --cpus-per-task="$THREADS" --cpu-bind=cores --exact bash -lc "$1"
}

appt_exec() {
  srun_exec apptainer exec --cleanenv "${APPT_ENV_ARGS[@]}" "${BIND_ARR[@]}" "$SIF" "$@"
}

container_has() {
  local tool="$1"
  appt_exec bash -c "command -v \"$tool\" >/dev/null 2>&1"
}

resolve_runtime_mode() {
  case "$K2_RUNTIME_MODE" in
    direct)
      echo "direct"
      ;;
    micromamba)
      echo "micromamba"
      ;;
    auto)
      if container_has kraken2 && container_has bracken; then
        echo "direct"
      else
        echo "micromamba"
      fi
      ;;
    *)
      echo "[error] Invalid K2_RUNTIME_MODE=$K2_RUNTIME_MODE (use: auto|direct|micromamba)" >&2
      exit 1
      ;;
  esac
}

run_tool() {
  if [[ "$RUNTIME_MODE" == "direct" ]]; then
    appt_exec "$@"
  else
    appt_exec micromamba run -n "$MICROMAMBA_ENV_NAME" "$@"
  fi
}

# -----------------------------
# FASTQ detection
# -----------------------------
shopt -s nullglob

detect_fastq() {
  local pattern
  local suf
  local files=()

  for pattern in "$@"; do
    for suf in ".fastq.gz" ".fq.gz" ".fastq" ".fq"; do
      files=( "$DATA_DIR"/"${SAMPLE}"*"${pattern}"*"${suf}" )
      if ((${#files[@]} > 0)); then
        printf '%s\n' "${files[0]}"
        return 0
      fi
    done
  done
  return 1
}

detect_any_fastq() {
  local suf
  local files=()
  for suf in ".fastq.gz" ".fq.gz" ".fastq" ".fq"; do
    files=( "$DATA_DIR"/"${SAMPLE}"*"${suf}" )
    if ((${#files[@]} > 0)); then
      printf '%s\n' "${files[0]}"
      return 0
    fi
  done
  return 1
}

R1="${R1:-}"
R2="${R2:-}"

if [[ -z "$R1" ]]; then
  R1="$(detect_fastq "_1" "_R1" ".R1" "R1" || detect_any_fastq || true)"
fi

[[ -f "$R1" ]] || {
  echo "[error] R1 not found. DATA_DIR=$DATA_DIR SAMPLE=$SAMPLE" >&2
  exit 1
}

if [[ "$PAIRED" == "true" ]]; then
  if [[ -z "$R2" ]]; then
    R2="$(detect_fastq "_2" "_R2" ".R2" "R2" || true)"
  fi
  [[ -f "${R2:-}" ]] || {
    echo "[error] PAIRED=true but R2 not found. DATA_DIR=$DATA_DIR SAMPLE=$SAMPLE" >&2
    exit 1
  }
fi

shopt -u nullglob

# -----------------------------
# Runtime mode
# -----------------------------
RUNTIME_MODE="$(resolve_runtime_mode)"

# If micromamba is required, isolate cache/lock per job
if [[ "$RUNTIME_MODE" == "micromamba" ]]; then
  XDG_CACHE_INNER="/tmp/xdg_cache_${USER}_${SLURM_JOB_ID:-$$}"
  APPT_ENV_ARGS+=(--env "XDG_CACHE_HOME=$XDG_CACHE_INNER")
  APPT_ENV_ARGS+=(--env "MAMBA_NO_BANNER=1")
fi

# -----------------------------
# Logging
# -----------------------------
echo "[info] DATASET=$DATASET SAMPLE=$SAMPLE"
echo "[info] R1=$R1"
if [[ "$PAIRED" == "true" ]]; then
  echo "[info] R2=$R2 (paired)"
else
  echo "[info] single-end"
fi
echo "[info] OUT_DIR=$OUT_DIR"
echo "[info] THREADS=$THREADS MEM_GB=$MEM_GB READ_LEN=$READ_LEN BRACKEN_THRESHOLD=$BRACKEN_THRESHOLD"
echo "[info] RUNTIME_MODE=$RUNTIME_MODE"

{
  echo "=== Runtime mode ==="
  echo "$RUNTIME_MODE"
  echo "=== Kraken2 version ==="
  run_tool kraken2 --version 2>&1 | sed -n '1p' || true
  echo "=== Bracken version ==="
  run_tool bracken -v 2>&1 | sed -n '1p' || true
  echo "=== Timestamp ==="
  date -Iseconds
} > "$VERSIONS" 2>&1 || true

cat > "$PARAMS_JSON" <<EOF
{"sample":"$SAMPLE","dataset":"$DATASET","paired":"$PAIRED",
 "r1":"$R1","r2":"${R2:-}","threads":$THREADS,"mem_gb":$MEM_GB,"bracken_threshold":$BRACKEN_THRESHOLD,
 "kraken_db":"$K2_DB","read_len":$READ_LEN,"runtime_mode":"$RUNTIME_MODE"}
EOF

# -----------------------------
# Optional QC read counting
# -----------------------------
if [[ "$KEEP_QC" == "1" ]]; then
  if [[ "$PAIRED" == "true" ]]; then
    L1=$(srun_shell 'if command -v pigz >/dev/null 2>&1; then pigz -dc "'"$R1"'"; else zcat -f "'"$R1"'"; fi | wc -l')
    T1=$((L1 / 4))

    L2=$(srun_shell 'if command -v pigz >/dev/null 2>&1; then pigz -dc "'"$R2"'"; else zcat -f "'"$R2"'"; fi | wc -l')
    T2=$((L2 / 4))

    TOTAL=$((T1 + T2))
  else
    L1=$(srun_shell 'if command -v pigz >/dev/null 2>&1; then pigz -dc "'"$R1"'"; else zcat -f "'"$R1"'"; fi | wc -l')
    TOTAL=$((L1 / 4))
  fi
  echo "$TOTAL" > "$TOTAL_FILE"
fi

# -----------------------------
# Kraken2
# -----------------------------
if [[ "$PAIRED" == "true" ]]; then
  run_tool kraken2 \
    --db "$K2_DB" \
    --threads "$THREADS" \
    --paired "$R1" "$R2" \
    --report "$KRAKEN_REPORT" \
    --output "$KRAKEN_OUT"
else
  run_tool kraken2 \
    --db "$K2_DB" \
    --threads "$THREADS" \
    "$R1" \
    --report "$KRAKEN_REPORT" \
    --output "$KRAKEN_OUT"
fi

# -----------------------------
# Bracken species
# -----------------------------
run_tool bracken \
  -d "$K2_DB" \
  -i "$KRAKEN_REPORT" \
  -o "$BRACKEN_S" \
  -r "$READ_LEN" \
  -l S \
  -t "$BRACKEN_THRESHOLD"

# -----------------------------
# Optional Bracken genus
# -----------------------------
if [[ "$KEEP_BRACKEN_GENUS" == "1" ]]; then
  run_tool bracken \
    -d "$K2_DB" \
    -i "$KRAKEN_REPORT" \
    -o "$BRACKEN_G" \
    -r "$READ_LEN" \
    -l G \
    -t "$BRACKEN_THRESHOLD"
fi

# -----------------------------
# Cleanup
# -----------------------------
if [[ "$KEEP_KRAKEN_REPORT" != "1" ]]; then
  rm -f "$KRAKEN_REPORT"
fi

if [[ "$KEEP_BRACKEN_GENUS" != "1" ]]; then
  rm -f "$BRACKEN_G"
fi

if [[ "$KEEP_KRAKEN_OUT" != "1" ]]; then
  rm -f "$OUT_DIR/${SAMPLE}.kraken2.out" 2>/dev/null || true
fi

echo "[ok] Kraken2+Bracken complete for $SAMPLE"
