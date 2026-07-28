#!/bin/bash
#SBATCH -J metaphlan4_profile
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o metaphlan4.%x.%j.out
#SBATCH -e metaphlan4.%x.%j.err
#SBATCH --exclude=compute-[1,5,9,13,17,21]

set -euo pipefail
IFS=$'\n\t'
umask 0022
export LC_ALL=C

# -----------------------------
# Resolve ROOT and load globals
# -----------------------------
if [[ -z "${ROOT:-}" || ! -d "$ROOT" ]]; then
  ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
fi
source "$ROOT/config/global.env" || true

# -----------------------------
# Required exported vars
# -----------------------------
: "${SAMPLE:?SAMPLE not set}"
: "${DATASET:?DATASET not set}"
: "${DATA_DIR:?DATA_DIR not set}"
PAIRED="${PAIRED:-false}"

: "${SIF:?SIF not set}"
: "${BIND:?BIND not set}"
: "${MPA_DB:?MPA_DB not set}"

# -----------------------------
# Tunables
# -----------------------------
THREADS="${MPA_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

# Keep the MetaPhlAn intermediate mapping file?
# Supports old env name KEEP_BOWTIE2OUT for compatibility.
KEEP_MAPOUT="${KEEP_MAPOUT:-${KEEP_BOWTIE2OUT:-0}}"

# Count total reads for QC/open-world downstream?
KEEP_QC="${KEEP_QC:-0}"

OUT_ROOT="${OUT_ROOT:-$ROOT/results}"
OUT_DIR="$OUT_ROOT/${DATASET}/metaphlan4"
mkdir -p "$OUT_DIR"

RAW="$OUT_DIR/${SAMPLE}.metaphlan.tsv"
TOTAL_FILE="$OUT_DIR/${SAMPLE}.total_reads.txt"
VERSIONS="$OUT_DIR/${SAMPLE}.versions.txt"
PARAMS_JSON="$OUT_DIR/${SAMPLE}.params.json"

# -----------------------------
# Per-job isolated HOME/cache
# Avoid micromamba/home lock contention
# -----------------------------
JOB_TMP_BASE="${SLURM_TMPDIR:-/tmp}/${USER:-user}.metaphlan.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
JOB_HOME="$JOB_TMP_BASE/home"
mkdir -p "$JOB_HOME/.cache" "$JOB_HOME/.config" "$JOB_HOME/.local"

cleanup() {
  # Remove temporary mapout if we created one
  if [[ -n "${TMP_MAPOUT_TO_DELETE:-}" && -f "${TMP_MAPOUT_TO_DELETE:-}" ]]; then
    rm -f "$TMP_MAPOUT_TO_DELETE" || true
  fi
  # Remove temp home
  rm -rf "$JOB_TMP_BASE" || true
}
trap cleanup EXIT

# -----------------------------
# Build apptainer command
# -----------------------------
APPTAINER_CMD=(apptainer exec --cleanenv --home "$JOB_HOME")
if [[ -n "${BIND:-}" ]]; then
  APPTAINER_CMD+=(-B "$BIND")
fi
APPTAINER_CMD+=("$SIF")

# -----------------------------
# FASTQ autodetect
# -----------------------------
detect_fastq() {
  local pattern="$1"
  local f
  shopt -s nullglob
  for suf in ".fastq.gz" ".fq.gz" ".fastq" ".fq"; do
    for f in "$DATA_DIR"/"${SAMPLE}"*"$pattern"*"$suf"; do
      [[ -e "$f" ]] && {
        printf '%s\n' "$f"
        shopt -u nullglob
        return 0
      }
    done
  done
  shopt -u nullglob
  return 1
}

R1="${R1:-}"
R2="${R2:-}"

[[ -z "$R1" ]] && R1="$(detect_fastq "_1" || detect_fastq "" || true)"
[[ -f "$R1" ]] || {
  echo "ERROR: R1 not found. DATA_DIR=$DATA_DIR SAMPLE=$SAMPLE" >&2
  exit 1
}

if [[ "$PAIRED" == "true" ]]; then
  [[ -z "$R2" ]] && R2="$(detect_fastq "_2" || true)"
  [[ -f "$R2" ]] || {
    echo "ERROR: PAIRED=true but R2 not found. DATA_DIR=$DATA_DIR SAMPLE=$SAMPLE" >&2
    exit 1
  }
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
echo "[info] OUT_DIR=$OUT_DIR THREADS=$THREADS"
echo "[info] KEEP_MAPOUT=$KEEP_MAPOUT KEEP_QC=$KEEP_QC"
echo "[info] JOB_HOME=$JOB_HOME"
echo "[info] Host node: $(hostname)"

# -----------------------------
# Versions
# -----------------------------
{
  echo "=== MetaPhlAn version ==="
  "${APPTAINER_CMD[@]}" bash -lc "
    export HOME='$JOB_HOME'
    export XDG_CACHE_HOME='$JOB_HOME/.cache'
    export XDG_CONFIG_HOME='$JOB_HOME/.config'
    export XDG_DATA_HOME='$JOB_HOME/.local/share'
    micromamba run -n taxonomic_tools metaphlan --version 2>&1 | head -n 2
  " || true
  echo "=== Timestamp ==="
  date -Iseconds
} > "$VERSIONS" 2>&1 || true

# -----------------------------
# Optional total-read counting
# -----------------------------
if [[ "$KEEP_QC" == "1" ]]; then
  if [[ "$PAIRED" == "true" ]]; then
    if command -v pigz >/dev/null 2>&1; then
      L1=$(pigz -dc "$R1" | wc -l)
      L2=$(pigz -dc "$R2" | wc -l)
    else
      L1=$(zcat -f "$R1" | wc -l)
      L2=$(zcat -f "$R2" | wc -l)
    fi
    T1=$((L1 / 4))
    T2=$((L2 / 4))
    TOTAL=$((T1 + T2))
  else
    if command -v pigz >/dev/null 2>&1; then
      L1=$(pigz -dc "$R1" | wc -l)
    else
      L1=$(zcat -f "$R1" | wc -l)
    fi
    TOTAL=$((L1 / 4))
  fi

  echo "$TOTAL" > "$TOTAL_FILE"
  echo "[info] FASTQ total reads: $TOTAL"
else
  rm -f "$TOTAL_FILE"
fi

# -----------------------------
# Input / mapout handling
# MetaPhlAn expects --mapout for paired comma-separated FASTQs
# -----------------------------
if [[ "$PAIRED" == "true" ]]; then
  MPA_INPUT="${R1},${R2}"

  if [[ "$KEEP_MAPOUT" == "1" ]]; then
    MAPOUT="$OUT_DIR/${SAMPLE}.mapout.bz2"
  else
    MAPOUT="${SLURM_TMPDIR:-/tmp}/${SAMPLE}.mapout.bz2"
    TMP_MAPOUT_TO_DELETE="$MAPOUT"
  fi
else
  MPA_INPUT="$R1"
  MAPOUT=""
fi

# -----------------------------
# Params record
# -----------------------------
cat > "$PARAMS_JSON" <<EOF
{"sample":"$SAMPLE","dataset":"$DATASET","paired":"$PAIRED",
 "r1":"$R1","r2":"${R2:-}","threads":$THREADS,
 "metaphlan_db":"$MPA_DB","keep_mapout":"$KEEP_MAPOUT","keep_qc":"$KEEP_QC",
 "output_tsv":"$RAW","mapout":"${MAPOUT:-}"}
EOF

# -----------------------------
# Run MetaPhlAn
# -----------------------------
RUN_SCRIPT=$(cat <<EOF
set -euo pipefail
export HOME='$JOB_HOME'
export XDG_CACHE_HOME='$JOB_HOME/.cache'
export XDG_CONFIG_HOME='$JOB_HOME/.config'
export XDG_DATA_HOME='$JOB_HOME/.local/share'
mkdir -p "\$HOME/.cache" "\$HOME/.config" "\$HOME/.local/share"

cmd=(
  micromamba run -n taxonomic_tools
  metaphlan
  "$MPA_INPUT"
  --input_type fastq
  --db_dir "$MPA_DB"
  --nproc "$THREADS"
  --ignore_eukaryotes
  --ignore_archaea
)

if [[ "$PAIRED" == "true" ]]; then
  cmd+=( --mapout "$MAPOUT" )
fi

cmd+=( -o "$RAW" )

echo "[cmd] \${cmd[*]}" >&2
"\${cmd[@]}"
EOF
)

echo "[info] Launching MetaPhlAn..."
"${APPTAINER_CMD[@]}" bash -lc "$RUN_SCRIPT"

[[ -s "$RAW" ]] || {
  echo "ERROR: MetaPhlAn finished but output file is missing/empty: $RAW" >&2
  exit 1
}

# -----------------------------
# Done
# -----------------------------
echo "[done] MetaPhlAn output: $RAW"
if [[ "$PAIRED" == "true" && "$KEEP_MAPOUT" == "1" ]]; then
  echo "[done] MetaPhlAn mapout: $MAPOUT"
fi
echo "[ok] MetaPhlAn complete for $SAMPLE"