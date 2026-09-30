#!/bin/bash
#SBATCH --job-name=bracken_denominator
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=bracken_denominator.%j.out
#SBATCH --error=bracken_denominator.%j.err
set -euo pipefail
: "${PROJECT:?Export the repository root}"
: "${AUDIT_ROOT:?Export a fresh output directory}"
cd "$PROJECT"
for spec in yachida:yachida_strict_final_20260823 feng:feng_strict_final_20260911 zeller:zeller_strict_final_20260911; do
    IFS=: read -r cohort run_name <<< "$spec"
    python3 -B analysis_v2/scripts/audit_bracken_denominators.py \
        --cohort "$cohort" \
        --state-dir "$PROJECT/work/$run_name/state" \
        --results-root "$PROJECT/work/$run_name/results" \
        --outdir "$AUDIT_ROOT/$cohort"
done
