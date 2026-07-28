#!/usr/bin/env bash
#SBATCH -J make_spike_pool
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o %x.%A.out
#SBATCH -e %x.%A.err

set -euo pipefail
IFS=$'\n\t'

: "${IMG:?set IMG (taxonomic-tools.sif)}"
: "${REF_FASTA:?set REF_FASTA (reference genome fasta)}"
: "${WORK:?set WORK (output dir)}"

LABEL="${LABEL:-Bfrag}"
PAIRED="${PAIRED:-true}"
READLEN="${READLEN:-150}"
FRAGMEAN="${FRAGMEAN:-350}"
FRAGSD="${FRAGSD:-10}"
COVERAGE="${COVERAGE:-2000}"
MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"

mkdir -p "$WORK/logs"
cd "$WORK"

# Unique run identifier for temporary ART outputs. This prevents failed/interrupted
# jobs from leaving final-looking files that confuse later reruns.
RUN_ID="${SLURM_JOB_ID:-manual_$$}"

# Run container work. Arguments are passed positionally to avoid fragile quoting
# inside bash -lc strings.
singularity exec \
  -B "$WORK:$WORK" \
  -B "$(dirname "$REF_FASTA")":"$(dirname "$REF_FASTA")" \
  "$IMG" bash -s -- \
    "$WORK" "$REF_FASTA" "$LABEL" "$PAIRED" "$READLEN" "$FRAGMEAN" "$FRAGSD" "$COVERAGE" "$MAMBA_ENV" "$RUN_ID" <<'EOS'
set -euo pipefail
IFS=$'\n\t'

WORK="$1"
REF_FASTA="$2"
LABEL="$3"
PAIRED="$4"
READLEN="$5"
FRAGMEAN="$6"
FRAGSD="$7"
COVERAGE="$8"
MAMBA_ENV="$9"
RUN_ID="${10}"

cd "$WORK"

if [[ "$PAIRED" == "true" ]]; then
  raw1="${LABEL}.raw_1.fq"
  raw2="${LABEL}.raw_2.fq"
  tmp_prefix="${LABEL}.art.${RUN_ID}"
  tmp1="${tmp_prefix}1.fq"
  tmp2="${tmp_prefix}2.fq"

  if [[ -s "$raw1" && -s "$raw2" ]]; then
    echo "[SKIP] Raw paired files already exist for ${LABEL}: ${raw1} ${raw2}"
  elif [[ -e "$raw1" || -e "$raw2" ]]; then
    echo "[ERROR] Incomplete raw paired files for ${LABEL}. Remove both ${raw1}/${raw2} and rerun." >&2
    exit 1
  else
    rm -f "$tmp1" "$tmp2"

    micromamba run -n "$MAMBA_ENV" art_illumina \
      -ss HS25 -i "$REF_FASTA" -l "$READLEN" -f "$COVERAGE" -p \
      -m "$FRAGMEAN" -s "$FRAGSD" -na -o "$tmp_prefix"

    test -s "$tmp1"
    test -s "$tmp2"

    mv -f "$tmp1" "$raw1"
    mv -f "$tmp2" "$raw2"
    echo "[OK] Created raw paired files for ${LABEL}: ${raw1} ${raw2}"
  fi
else
  raw="${LABEL}.raw.fq"
  tmp_prefix="${LABEL}.art.${RUN_ID}"
  tmp="${tmp_prefix}.fq"

  if [[ -s "$raw" ]]; then
    echo "[SKIP] Raw single-end file already exists for ${LABEL}: ${raw}"
  elif [[ -e "$raw" ]]; then
    echo "[ERROR] Incomplete raw single-end file for ${LABEL}. Remove ${raw} and rerun." >&2
    exit 1
  else
    rm -f "$tmp"

    micromamba run -n "$MAMBA_ENV" art_illumina \
      -ss HS25 -i "$REF_FASTA" -l "$READLEN" -f "$COVERAGE" -na -o "$tmp_prefix"

    test -s "$tmp"

    mv -f "$tmp" "$raw"
    echo "[OK] Created raw single-end file for ${LABEL}: ${raw}"
  fi
fi
EOS

# If sbatch -o/-e was not overridden by spikein_prepare_pools.sh, SLURM writes
# logs in the submit directory. Move those default logs into WORK/logs when they exist.
OUTLOG="$(pwd)/logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"
ERRLOG="$(pwd)/logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err"
SUBMIT_OUT="${SLURM_SUBMIT_DIR:-$(pwd)}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out"
SUBMIT_ERR="${SLURM_SUBMIT_DIR:-$(pwd)}/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err"

[[ -f "$SUBMIT_OUT" && "$SUBMIT_OUT" != "$OUTLOG" ]] && mv -f "$SUBMIT_OUT" "$OUTLOG" || true
[[ -f "$SUBMIT_ERR" && "$SUBMIT_ERR" != "$ERRLOG" ]] && mv -f "$SUBMIT_ERR" "$ERRLOG" || true

echo "[OK] Spike pool simulation finished in $WORK"
echo "[OK] Default logs, if not overridden by sbatch -o/-e: $OUTLOG  $ERRLOG"