#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage(){
  cat <<'USAGE'
Usage:
  spikein_prepare_pools.sh --env spikein.env

Submits, for each taxon in SPIKE_PANEL:
  1) ART pool simulation (make_spike_pool.sbatch)
  2) fastp preprocessing (fastp_spike_pool.sbatch) dependent on (1)

Assumes spike_panel.tsv has columns:
  label  taxon_name  assembly  fasta  weight  url

This script does NOT download FASTAs.
It expects `fasta` paths in SPIKE_PANEL to already exist on disk.

Writes:
  ${POOLS_DIR:-$WORK/pools}/jobs.tsv
  ${POOLS_DIR:-$WORK/pools}/pool_settings.tsv
  ${POOLS_DIR:-$WORK/pools}/logs/<label>.make.<jobid>.{out,err}
  ${POOLS_DIR:-$WORK/pools}/logs/<label>.fastp.<jobid>.{out,err}

Important behavior:
  - If <LABEL>.pool_1.fq and <LABEL>.pool_2.fq already exist and are non-empty,
    this script skips that label.
  - If only one mate exists, the script stops with an error.
  - If pool files exist but are older than raw files, the script stops with an error,
    because the directory may contain mixed/stale outputs from multiple runs.
  - Use a clean POOLS_DIR for each parameter set, e.g. pools_readlen100_cov2000.
USAGE
}

ENV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$ENV" ]] || { usage; exit 2; }
[[ -f "$ENV" ]] || { echo "[ERROR] Missing env file: $ENV" >&2; exit 1; }

# shellcheck disable=SC1090
source "$ENV"

: "${IMG:?set IMG in env}"
: "${WORK:?set WORK in env}"
: "${SPIKE_PANEL:?set SPIKE_PANEL in env}"

command -v sbatch >/dev/null 2>&1 || {
  echo "[ERROR] sbatch not found. Run on a SLURM login node." >&2
  exit 1
}

READLEN="${READLEN:-150}"
COVERAGE="${COVERAGE:-2000}"
FRAGMEAN="${FRAGMEAN:-350}"
FRAGSD="${FRAGSD:-10}"
THREADS="${THREADS:-4}"
MAMBA_ENV="${MAMBA_ENV:-taxonomic_tools}"

# Recommended: set this explicitly in spikein.env, for example:
# POOLS_DIR=$WORK/pools_readlen100_cov2000
POOLS_DIR="${POOLS_DIR:-$WORK/pools}"
mkdir -p "$POOLS_DIR/logs"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -s "$SCRIPT_DIR/make_spike_pool.sbatch" ]] || {
  echo "[ERROR] Missing make script: $SCRIPT_DIR/make_spike_pool.sbatch" >&2
  exit 1
}
[[ -s "$SCRIPT_DIR/fastp_spike_pool.sbatch" ]] || {
  echo "[ERROR] Missing fastp script: $SCRIPT_DIR/fastp_spike_pool.sbatch" >&2
  exit 1
}

JOBS_TSV="$POOLS_DIR/jobs.tsv"
SETTINGS_TSV="$POOLS_DIR/pool_settings.tsv"

echo -e "label\tmake_jobid\tfastp_jobid\tstatus" > "$JOBS_TSV"
cat > "$SETTINGS_TSV" <<SETTINGS
parameter	value
ENV	$ENV
IMG	$IMG
WORK	$WORK
POOLS_DIR	$POOLS_DIR
SPIKE_PANEL	$SPIKE_PANEL
READLEN	$READLEN
COVERAGE	$COVERAGE
FRAGMEAN	$FRAGMEAN
FRAGSD	$FRAGSD
THREADS	$THREADS
MAMBA_ENV	$MAMBA_ENV
SCRIPT_DIR	$SCRIPT_DIR
SETTINGS

check_pair_state(){
  local a="$1"
  local b="$2"
  local label="$3"
  local what="$4"

  if [[ -s "$a" && ! -s "$b" ]]; then
    echo "[ERROR] $label has only one $what mate: $a exists but $b is missing/empty." >&2
    exit 1
  fi
  if [[ ! -s "$a" && -s "$b" ]]; then
    echo "[ERROR] $label has only one $what mate: $b exists but $a is missing/empty." >&2
    exit 1
  fi
}

# Expect header:
# label taxon_name assembly fasta weight url
while IFS=$'\t' read -r label taxon_name assembly fasta weight url extra; do
  [[ -z "${label:-}" ]] && continue
  [[ "$label" =~ ^# ]] && continue

  if [[ -n "${extra:-}" ]]; then
    echo "[WARN] Extra columns found for label=$label in SPIKE_PANEL; ignoring columns after 'url'." >&2
  fi

  # FASTA must exist. No downloading here.
  if [[ -z "${fasta:-}" || "$fasta" == "-" ]]; then
    echo "[ERROR] Missing fasta path for label=$label in SPIKE_PANEL column 'fasta'." >&2
    exit 1
  fi
  [[ -s "$fasta" ]] || {
    echo "[ERROR] FASTA missing/empty for $label: $fasta" >&2
    exit 1
  }

  raw1="$POOLS_DIR/${label}.raw_1.fq"
  raw2="$POOLS_DIR/${label}.raw_2.fq"
  pool1="$POOLS_DIR/${label}.pool_1.fq"
  pool2="$POOLS_DIR/${label}.pool_2.fq"

  # Leftovers from the old ART prefix pattern usually mean an interrupted run.
  # Do not silently mix them with new outputs.
  if compgen -G "$POOLS_DIR/${label}.art*.fq" >/dev/null; then
    echo "[ERROR] Found leftover ART files for $label in $POOLS_DIR:" >&2
    ls -lh "$POOLS_DIR/${label}.art"*.fq >&2 || true
    echo "[ERROR] Move/delete these leftovers or use a clean POOLS_DIR before rerunning." >&2
    exit 1
  fi

  check_pair_state "$raw1" "$raw2" "$label" "raw"
  check_pair_state "$pool1" "$pool2" "$label" "pool"

  # If final pools already exist, skip, but protect against stale pools from a mixed rerun.
  if [[ -s "$pool1" && -s "$pool2" ]]; then
    if [[ -s "$raw1" && "$pool1" -ot "$raw1" ]]; then
      echo "[ERROR] $label pool files exist but are older than raw files." >&2
      echo "[ERROR] This suggests mixed/stale outputs in $POOLS_DIR." >&2
      echo "[ERROR] Use a clean POOLS_DIR or remove both raw and pool files for this label." >&2
      exit 1
    fi
    if [[ -s "$raw2" && "$pool2" -ot "$raw2" ]]; then
      echo "[ERROR] $label pool files exist but are older than raw files." >&2
      echo "[ERROR] This suggests mixed/stale outputs in $POOLS_DIR." >&2
      echo "[ERROR] Use a clean POOLS_DIR or remove both raw and pool files for this label." >&2
      exit 1
    fi

    echo -e "${label}\tSKIP\tSKIP\tpool_exists" >> "$JOBS_TSV"
    echo "[SKIP] $label pools already exist: $pool1 $pool2"
    continue
  fi

  make_jid=$(sbatch --parsable \
    --job-name="make_${label}" \
    --output="$POOLS_DIR/logs/${label}.make.%j.out" \
    --error="$POOLS_DIR/logs/${label}.make.%j.err" \
    --export=ALL,IMG="$IMG",WORK="$POOLS_DIR",REF_FASTA="$fasta",LABEL="$label",PAIRED=true,READLEN="$READLEN",COVERAGE="$COVERAGE",FRAGMEAN="$FRAGMEAN",FRAGSD="$FRAGSD",MAMBA_ENV="$MAMBA_ENV" \
    "$SCRIPT_DIR/make_spike_pool.sbatch")

  fastp_jid=$(sbatch --parsable \
    --job-name="fastp_${label}" \
    --output="$POOLS_DIR/logs/${label}.fastp.%j.out" \
    --error="$POOLS_DIR/logs/${label}.fastp.%j.err" \
    --dependency=afterok:"$make_jid" \
    --export=ALL,IMG="$IMG",WORK="$POOLS_DIR",LABEL="$label",THREADS="$THREADS",MAMBA_ENV="$MAMBA_ENV" \
    "$SCRIPT_DIR/fastp_spike_pool.sbatch")

  echo -e "${label}\t${make_jid}\t${fastp_jid}\tsubmitted" >> "$JOBS_TSV"
  echo "[SUBMITTED] $label make=$make_jid fastp=$fastp_jid"
done < <(tail -n +2 "$SPIKE_PANEL")

echo "[OK] Wrote $JOBS_TSV"
echo "[OK] Wrote $SETTINGS_TSV"
echo "[OK] Logs will be written to $POOLS_DIR/logs"