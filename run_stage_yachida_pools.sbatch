#!/usr/bin/env bash
#SBATCH --job-name=stage_yachida_pools
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=stage_yachida_pools.%j.out
#SBATCH --error=stage_yachida_pools.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${SOURCE_POOLS_DIR:?export SOURCE_POOLS_DIR}"
: "${DESTINATION_POOLS_DIR:?export DESTINATION_POOLS_DIR}"
bash "$ROOT/spikes/scripts/spikein/stage_finalized_pools.sh" \
  --source "$SOURCE_POOLS_DIR" \
  --destination "$DESTINATION_POOLS_DIR"
