#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the MaAsLin2 analysis container}"

echo "[1/2] Upstream spike-generation and profiling preflight"
GLOBAL_ENV="${GLOBAL_ENV:-$ROOT/config/global.env}" \
SPIKE_ENV="${SPIKE_ENV:-$ROOT/spikes/spikein.env}" \
bash scripts/preflight.sh

echo "[2/2] Downstream statistical-analysis preflight"
PROJECT="$ROOT" \
SIF="$ANALYSIS_SIF" \
bash cluster_preflight_original_unpaired.sh

echo "[PASS] Complete repository preflight passed."
