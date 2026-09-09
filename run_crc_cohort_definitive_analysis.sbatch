#!/usr/bin/env bash
#SBATCH --job-name=crc_v2_final
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=crc_v2_final.%j.out
#SBATCH --error=crc_v2_final.%j.err
set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"; cd "$ROOT"
bash analysis_v2/run_crc_cohort_definitive_analysis.sh
