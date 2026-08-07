#!/usr/bin/env bash
#SBATCH --job-name=index_spike_pools
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=index_spike_pools.%j.out
#SBATCH --error=index_spike_pools.%j.err

set -euo pipefail
ROOT="${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
cd "$ROOT"
: "${SPIKE_ENV:?export SPIKE_ENV}"
bash "$ROOT/spikes/scripts/spikein/index_finalized_pools.sh" --env "$SPIKE_ENV"
