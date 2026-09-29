#!/usr/bin/env bash
#SBATCH --job-name=cohort_v2_final
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=cohort_v2_final.%j.out
#SBATCH --error=cohort_v2_final.%j.err
set -euo pipefail
ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/analysis_v2/run_cohort_definitive_analysis.sh" ]] || {
  echo "[ERROR] Submit from the project root or export PROJECT" >&2
  exit 1
}
cd "$ROOT"
bash analysis_v2/run_cohort_definitive_analysis.sh
