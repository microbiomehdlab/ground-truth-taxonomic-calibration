#!/usr/bin/env bash
#SBATCH --job-name=legacy_biomarker_dev
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

set -euo pipefail

ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/analysis_v2/run_legacy_biomarker_development.sh" ]] || {
  echo "[ERROR] Submit this job from the project root or export PROJECT" >&2
  exit 1
}
cd "$ROOT"

: "${DEV_INPUT_ROOT:?Export DEV_INPUT_ROOT before submitting}"
: "${ANALYSIS_SIF:?Export ANALYSIS_SIF before submitting}"
: "${OUTDIR:?Export OUTDIR before submitting}"

bash analysis_v2/run_legacy_biomarker_development.sh
