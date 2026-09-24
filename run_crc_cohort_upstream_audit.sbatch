#!/usr/bin/env bash
#SBATCH --job-name=crc_upstream_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=crc_upstream_audit.%j.out
#SBATCH --error=crc_upstream_audit.%j.err

set -euo pipefail
IFS=$'\n\t'

ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/analysis_v2/scripts/seal_crc_cohort_upstream.py" ]] || {
  echo "[ERROR] Submit from the project root or export PROJECT" >&2
  exit 1
}
cd "$ROOT"

: "${CRC_ENV:?export CRC_ENV}"
: "${COHORT:?set COHORT to feng or zeller}"
: "${PRODUCTION_MANIFEST:?set the frozen production manifest}"
: "${INDEPENDENT_MANIFEST:?set the frozen independent-subset manifest}"
: "${EXPECTED_SAMPLES:?set the frozen cohort size}"

[[ "$COHORT" == feng || "$COHORT" == zeller ]] || {
  echo "[ERROR] COHORT must be feng or zeller" >&2
  exit 2
}
[[ "$EXPECTED_SAMPLES" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] EXPECTED_SAMPLES must be a positive integer" >&2
  exit 2
}

source "$CRC_ENV"
: "${CRC_STATE_DIR:?set in CRC_ENV}"
: "${CRC_SCRATCH_ROOT:?set in CRC_ENV}"
: "${PERSISTENT_RESULTS_ROOT:?set in CRC_ENV}"
: "${PERSISTENT_QC_ROOT:?set in CRC_ENV}"

python3 "$ROOT/analysis_v2/scripts/seal_crc_cohort_upstream.py" \
  --cohort "$COHORT" \
  --manifest "$PRODUCTION_MANIFEST" \
  --independent-manifest "$INDEPENDENT_MANIFEST" \
  --state-dir "$CRC_STATE_DIR" \
  --results-root "$PERSISTENT_RESULTS_ROOT" \
  --qc-root "$PERSISTENT_QC_ROOT" \
  --scratch-root "$CRC_SCRATCH_ROOT" \
  --expected-samples "$EXPECTED_SAMPLES" \
  --outdir "$CRC_STATE_DIR/production_seal"
