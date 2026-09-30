#!/bin/bash
#SBATCH --job-name=bracken_recovery
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=bracken_recovery.%j.out
#SBATCH --error=bracken_recovery.%j.err
set -euo pipefail
: "${PROJECT:?Export repository root}"
: "${AUDIT_ROOT:?Export successful denominator audit root}"
: "${RECOVERY_ROOT:?Export fresh recovery output root}"
: "${PREFLIGHT_STAMP:?Export timestamp of the shared preflight runs}"
cd "$PROJECT"
for cohort in yachida feng zeller; do
    canonical="$PROJECT/work/analysis_v2_${cohort}_shared_preflight_${PREFLIGHT_STAMP}/canonical/canonical_input.tsv"
    test -s "$canonical" || { echo "Missing canonical input: $canonical" >&2; exit 1; }
    python3 -B analysis_v2/scripts/derive_bracken_all_input_recovery.py \
        --cohort "$cohort" --project "$PROJECT" \
        --audit "$AUDIT_ROOT/$cohort" --canonical "$canonical" \
        --community-allocation frozen-panel-reconstruction \
        --outdir "$RECOVERY_ROOT/$cohort"
done
