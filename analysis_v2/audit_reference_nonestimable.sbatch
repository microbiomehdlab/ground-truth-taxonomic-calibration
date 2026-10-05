#!/usr/bin/env bash
#SBATCH --job-name=reference_ne_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=reference_ne_audit.%j.out
#SBATCH --error=reference_ne_audit.%j.err
set -euo pipefail
: "${PROJECT:?}" "${REFERENCE_ROOT:?}" "${REFERENCE_NE_ROOT:?}"
cd "$PROJECT"
python3 -B analysis_v2/scripts/audit_reference_nonestimable.py \
    --root "$REFERENCE_ROOT" --report "$REFERENCE_ROOT/REPORT" --out "$REFERENCE_NE_ROOT"
