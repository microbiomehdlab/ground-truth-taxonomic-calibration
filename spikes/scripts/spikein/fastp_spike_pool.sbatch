#!/usr/bin/env bash
#SBATCH -J fastp_spike_pool
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o %x.%A.out
#SBATCH -e %x.%A.err

set -euo pipefail
IFS=$'\n\t'

: "${IMG:?set IMG (taxonomic-tools.sif)}"
: "${WORK:?set WORK (dir containing <LABEL>.raw_1.fq and <LABEL>.raw_2.fq)}"

LABEL="${LABEL:-Bfrag}"
THREADS="${THREADS:-4}"
MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"

mkdir -p "$WORK/logs"
cd "$WORK"

RAW1="${LABEL}.raw_1.fq"
RAW2="${LABEL}.raw_2.fq"
POOL1="${LABEL}.pool_1.fq"
POOL2="${LABEL}.pool_2.fq"
HTML="logs/${LABEL}.fastp.html"
JSON="logs/${LABEL}.fastp.json"

RUN_ID="${SLURM_JOB_ID:-manual_$$}"
TMP1="${LABEL}.pool_1.tmp.${RUN_ID}.fq"
TMP2="${LABEL}.pool_2.tmp.${RUN_ID}.fq"
TMPHTML="logs/${LABEL}.fastp.tmp.${RUN_ID}.html"
TMPJSON="logs/${LABEL}.fastp.tmp.${RUN_ID}.json"

count_fastq_reads() {
  awk 'END{printf "%.0f", NR/4}' "$1"
}

# Raw paired files must both exist and contain the same number of records.
[[ -s "$RAW1" ]] || { echo "[ERROR] Missing/empty raw R1 for $LABEL: $WORK/$RAW1" >&2; exit 1; }
[[ -s "$RAW2" ]] || { echo "[ERROR] Missing/empty raw R2 for $LABEL: $WORK/$RAW2" >&2; exit 1; }

raw_n1=$(count_fastq_reads "$RAW1")
raw_n2=$(count_fastq_reads "$RAW2")
if [[ "$raw_n1" != "$raw_n2" ]]; then
  echo "[ERROR] Raw mate count mismatch for $LABEL: $RAW1=$raw_n1 $RAW2=$raw_n2" >&2
  exit 1
fi

# If final pools already exist, validate them and skip. Do not silently accept
# incomplete, stale, or unsynchronised paired-end pools.
if [[ -s "$POOL1" || -s "$POOL2" ]]; then
  if [[ ! -s "$POOL1" || ! -s "$POOL2" ]]; then
    echo "[ERROR] Incomplete existing pool for $LABEL. Remove both $POOL1 and $POOL2, or use a fresh POOLS_DIR." >&2
    exit 1
  fi

  if [[ "$POOL1" -ot "$RAW1" || "$POOL2" -ot "$RAW2" ]]; then
    echo "[ERROR] Existing pool for $LABEL is older than the raw FASTQ. This suggests a stale/mixed rerun." >&2
    echo "[ERROR] Use a fresh POOLS_DIR, or remove $POOL1 $POOL2 and rerun fastp for this label." >&2
    exit 1
  fi

  pool_n1=$(count_fastq_reads "$POOL1")
  pool_n2=$(count_fastq_reads "$POOL2")
  if [[ "$pool_n1" != "$pool_n2" ]]; then
    echo "[ERROR] Existing pool mate count mismatch for $LABEL: $POOL1=$pool_n1 $POOL2=$pool_n2" >&2
    echo "[ERROR] Use a fresh POOLS_DIR, or remove $POOL1 $POOL2 and rerun fastp for this label." >&2
    exit 1
  fi

  echo "[SKIP] Valid existing pool for $LABEL: $pool_n1 read pairs"
  exit 0
fi

# Clean only this job's temporary files, never final outputs.
rm -f "$TMP1" "$TMP2" "$TMPHTML" "$TMPJSON"

singularity exec -B "$WORK:$WORK" "$IMG" bash -lc "
  set -euo pipefail
  cd '$WORK'

  micromamba run -n '$MAMBA_ENV' fastp \
    -i '$RAW1' -I '$RAW2' \
    -o '$TMP1' -O '$TMP2' \
    --cut_front --cut_tail \
    --cut_window_size 4 --cut_mean_quality 20 \
    --trim_poly_g \
    --low_complexity_filter --complexity_threshold 30 \
    --length_required 60 \
    --detect_adapter_for_pe \
    --thread '$THREADS' \
    --html '$TMPHTML' \
    --json '$TMPJSON'
"

[[ -s "$TMP1" ]] || { echo "[ERROR] fastp produced missing/empty temporary R1 for $LABEL: $TMP1" >&2; exit 1; }
[[ -s "$TMP2" ]] || { echo "[ERROR] fastp produced missing/empty temporary R2 for $LABEL: $TMP2" >&2; exit 1; }

tmp_n1=$(count_fastq_reads "$TMP1")
tmp_n2=$(count_fastq_reads "$TMP2")
if [[ "$tmp_n1" != "$tmp_n2" ]]; then
  echo "[ERROR] fastp output mate count mismatch for $LABEL: $TMP1=$tmp_n1 $TMP2=$tmp_n2" >&2
  echo "[ERROR] Temporary files were left for inspection: $TMP1 $TMP2" >&2
  exit 1
fi

mv -f "$TMP1" "$POOL1"
mv -f "$TMP2" "$POOL2"
[[ -s "$TMPHTML" ]] && mv -f "$TMPHTML" "$HTML"
[[ -s "$TMPJSON" ]] && mv -f "$TMPJSON" "$JSON"

echo "[OK] Wrote $POOL1 and $POOL2 in $WORK"
echo "[OK] $LABEL pool contains $tmp_n1 synchronized read pairs"