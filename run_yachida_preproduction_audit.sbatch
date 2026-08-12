#!/usr/bin/env bash
#SBATCH --job-name=yachida_preflight
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=yachida_preflight.%j.out
#SBATCH --error=yachida_preflight.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
bash "$ROOT/run_yachida_preproduction_audit.sh"
