#!/usr/bin/env bash
#SBATCH --job-name=yachida_sampling_test
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=yachida_sampling_test.%j.out
#SBATCH --error=yachida_sampling_test.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${YACHIDA_ENV:?export YACHIDA_ENV}"
: "${SAMPLE_ID:?export SAMPLE_ID}"
export PROJECT="$ROOT"
bash "$ROOT/datasets/yachida/run_sampling_mode_integration_test.sh"
