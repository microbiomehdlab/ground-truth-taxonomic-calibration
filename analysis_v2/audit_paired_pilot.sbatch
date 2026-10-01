#!/bin/bash
#SBATCH --job-name=paired_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=paired_audit.%j.out
#SBATCH --error=paired_audit.%j.err
set -euo pipefail
: "${PROJECT:?}"
: "${ANALYSIS_SIF:?}"
: "${DA_PILOT_RESULTS:?}"
: "${PAIRED_AUDIT_ROOT:?}"
cd "$PROJECT"
test ! -e "$PAIRED_AUDIT_ROOT"
mkdir -p "$PAIRED_AUDIT_ROOT"
find "$DA_PILOT_RESULTS/attempts" -type f -print0 | sort -z | xargs -0 sha256sum > "$PAIRED_AUDIT_ROOT/input_files.sha256"
apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
    Rscript analysis_v2/scripts/audit_paired_pilot.R "$DA_PILOT_RESULTS" "$PAIRED_AUDIT_ROOT" "$PROJECT"
sha256sum -c --quiet "$PAIRED_AUDIT_ROOT/input_files.sha256"
git rev-parse HEAD > "$PAIRED_AUDIT_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$PAIRED_AUDIT_ROOT/analysis_image.sha256"
printf 'status\tPASS_DIAGNOSTIC_AUDIT\nmodels_refitted\t0\npolicy_changed\t0\n' > "$PAIRED_AUDIT_ROOT/SUCCESS"
(cd "$PAIRED_AUDIT_ROOT"; find . -maxdepth 1 -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Read-only paired diagnostics: $PAIRED_AUDIT_ROOT"
