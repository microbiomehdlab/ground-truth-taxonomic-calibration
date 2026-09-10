#!/usr/bin/env bash
#SBATCH --job-name=yachida_v2_final
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=yachida_v2_final.%j.out
#SBATCH --error=yachida_v2_final.%j.err
set -euo pipefail
ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/analysis_v2/run_yachida_definitive_analysis.sh" ]] || {
  echo "[ERROR] Submit from the project root or export PROJECT" >&2
  exit 1
}
cd "$ROOT"
bash analysis_v2/run_yachida_definitive_analysis.sh
