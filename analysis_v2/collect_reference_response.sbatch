#!/usr/bin/env bash
#SBATCH --job-name=reference_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=reference_collect.%j.out
#SBATCH --error=reference_collect.%j.err
set -euo pipefail
: "${REFERENCE_ROOT:?}"
cd "$REFERENCE_ROOT/source"
python3 -B analysis_v2/scripts/reference_response_pilot.py collect \
    --root "$REFERENCE_ROOT" --out "$REFERENCE_ROOT/REPORT"
