#!/usr/bin/env bash
#SBATCH --job-name=yachida_v2_final
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=yachida_v2_final.%j.out
#SBATCH --error=yachida_v2_final.%j.err
set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
bash analysis_v2/run_yachida_definitive_analysis.sh
