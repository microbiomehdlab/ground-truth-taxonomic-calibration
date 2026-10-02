#!/bin/bash
#SBATCH --job-name=da_validation_report
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=da_validation_report.%j.out
#SBATCH --error=da_validation_report.%j.err
set -euo pipefail
: "${PROJECT:?}" "${FINAL_VALIDATION_ROOT:?}"
cd "$PROJECT"
sha256sum -c --quiet "$FINAL_VALIDATION_ROOT/submission.sha256"
python3 -B analysis_v2/scripts/collect_final_validation.py "$FINAL_VALIDATION_ROOT"
