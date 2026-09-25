#!/usr/bin/env bash
#SBATCH --job-name=cohort_upstream_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=cohort_upstream_audit.%j.out
#SBATCH --error=cohort_upstream_audit.%j.err

# One public upstream audit for Yachida, Feng and Zeller. Every cohort is
# validated against the same contract and emits the same production_seal_v2
# layout; only boundary settings differ, and they are passed explicitly.
#
# The audit rehashes every retained output, so it must run in a compute job.
# It never recomputes a profile and never modifies the native source seal.
#
# Frozen production sizes and condition counts:
#   yachida  EXPECTED_SAMPLES=201  EXPECTED_CONDITIONS=Control=67,Adenoma=67,CRC=67
#   feng     EXPECTED_SAMPLES=154  EXPECTED_CONDITIONS=Control=61,Adenoma=47,CRC=46
#   zeller   EXPECTED_SAMPLES=156  EXPECTED_CONDITIONS=Control=61,Adenoma=42,CRC=53
# At those sizes the auditor cross-checks the counts against its frozen record
# and refuses a disagreement.

set -euo pipefail
IFS=$'\n\t'

ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/analysis_v2/scripts/seal_cohort_upstream.py" ]] || {
  echo "[ERROR] Submit from the project root or export PROJECT" >&2
  exit 1
}
cd "$ROOT"

: "${COHORT:?set COHORT to yachida, feng or zeller}"
: "${COHORT_ENV:?export COHORT_ENV to the cohort strict-production env file}"
: "${PRODUCTION_MANIFEST:?set the frozen production manifest}"
: "${INDEPENDENT_MANIFEST:?set the frozen independent-subset manifest}"
: "${SOURCE_SEAL:?set the verified native production seal to read provenance from}"
: "${EXPECTED_SAMPLES:?set the frozen cohort size}"
: "${EXPECTED_CONDITIONS:?set the frozen condition counts, e.g. Control=67,Adenoma=67,CRC=67}"
: "${OUTDIR:?set a new production_seal_v2 directory}"

case "$COHORT" in
  yachida|feng|zeller) ;;
  *) echo "[ERROR] COHORT must be yachida, feng or zeller" >&2; exit 2 ;;
esac
[[ "$EXPECTED_SAMPLES" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] EXPECTED_SAMPLES must be a positive integer" >&2
  exit 2
}

# shellcheck source=/dev/null
source "$COHORT_ENV"

# The state, scratch, results and QC roots are named differently per cohort
# environment; resolve them here so the auditor receives explicit paths.
STATE_DIR="${STATE_DIR:-${YACHIDA_STATE_DIR:-${CRC_STATE_DIR:-}}}"
SCRATCH_ROOT="${SCRATCH_ROOT:-${YACHIDA_SCRATCH_ROOT:-${CRC_SCRATCH_ROOT:-}}}"
: "${STATE_DIR:?set STATE_DIR or the cohort-specific state variable}"
: "${SCRATCH_ROOT:?set SCRATCH_ROOT or the cohort-specific scratch variable}"
: "${PERSISTENT_RESULTS_ROOT:?set in COHORT_ENV}"
: "${PERSISTENT_QC_ROOT:?set in COHORT_ENV}"

for required in "$PRODUCTION_MANIFEST" "$INDEPENDENT_MANIFEST"; do
  [[ -s "$required" ]] || { echo "[ERROR] missing or empty: $required" >&2; exit 2; }
done
[[ -s "$SOURCE_SEAL/SUCCESS" ]] || {
  echo "[ERROR] native source seal has no SUCCESS: $SOURCE_SEAL" >&2
  exit 2
}
[[ ! -e "$SOURCE_SEAL/AUDIT_IN_PROGRESS" ]] || {
  echo "[ERROR] native source seal is mid-audit: $SOURCE_SEAL" >&2
  exit 2
}

python3 "$ROOT/analysis_v2/scripts/seal_cohort_upstream.py" \
  --cohort "$COHORT" \
  --manifest "$PRODUCTION_MANIFEST" \
  --independent-manifest "$INDEPENDENT_MANIFEST" \
  --source-seal "$SOURCE_SEAL" \
  --state-dir "$STATE_DIR" \
  --results-root "$PERSISTENT_RESULTS_ROOT" \
  --qc-root "$PERSISTENT_QC_ROOT" \
  --scratch-root "$SCRATCH_ROOT" \
  --outdir "$OUTDIR" \
  --expected-samples "$EXPECTED_SAMPLES" \
  --expected-independent "${EXPECTED_INDEPENDENT:-30}" \
  --expected-conditions "$EXPECTED_CONDITIONS"

# The seal is only authoritative when SUCCESS, the checksum manifest and a
# complete verification are all present.
[[ -s "$OUTDIR/SUCCESS" && -s "$OUTDIR/production_seal.sha256" ]] || {
  echo "[ERROR] unified seal is incomplete: $OUTDIR" >&2
  exit 1
}
[[ ! -e "$OUTDIR/AUDIT_IN_PROGRESS" ]] || {
  echo "[ERROR] unified seal is still marked in progress: $OUTDIR" >&2
  exit 1
}
(cd "$OUTDIR" && sha256sum -c production_seal.sha256)
echo "[PASS] $COHORT unified upstream seal verified: $OUTDIR"
