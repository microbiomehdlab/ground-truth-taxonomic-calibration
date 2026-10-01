#!/bin/bash
#SBATCH --job-name=paired_refit_diag
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=paired_refit_diag.%j.out
#SBATCH --error=paired_refit_diag.%j.err
set -euo pipefail
: "${PROJECT:?}"
: "${ANALYSIS_SIF:?}"
: "${DA_PILOT_RESULTS:?}"
: "${PAIRED_AUDIT_ROOT:?}"
: "${PAIRED_REFIT_ROOT:?}"
cd "$PROJECT"
test ! -e "$PAIRED_REFIT_ROOT"
(cd "$PAIRED_AUDIT_ROOT"; sha256sum -c --quiet SHA256SUMS)
sha256sum -c --quiet "$PAIRED_AUDIT_ROOT/input_files.sha256"
mkdir -p "$PAIRED_REFIT_ROOT"
sha256sum "$PAIRED_AUDIT_ROOT/SHA256SUMS" "$PAIRED_AUDIT_ROOT/input_files.sha256" > "$PAIRED_REFIT_ROOT/input_provenance.sha256"
sha256sum analysis_v2/scripts/refit_paired_diagnostics.R analysis_v2/scripts/audit_paired_pilot.R \
    analysis_v2/refit_paired_diagnostics.sbatch > "$PAIRED_REFIT_ROOT/code.sha256"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
apptainer exec --cleanenv --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1 \
    --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF" \
    Rscript analysis_v2/scripts/refit_paired_diagnostics.R "$PAIRED_AUDIT_ROOT" "$DA_PILOT_RESULTS" "$PAIRED_REFIT_ROOT"
(cd "$PAIRED_AUDIT_ROOT"; sha256sum -c --quiet SHA256SUMS)
sha256sum -c --quiet "$PAIRED_AUDIT_ROOT/input_files.sha256"
sha256sum -c --quiet "$PAIRED_REFIT_ROOT/input_provenance.sha256"
sha256sum -c --quiet "$PAIRED_REFIT_ROOT/code.sha256"
git rev-parse HEAD > "$PAIRED_REFIT_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$PAIRED_REFIT_ROOT/analysis_image.sha256"
printf 'status\tPASS_DIAGNOSTIC_REFIT\nprimary_results_replaced\t0\npolicy_changed\t0\n' > "$PAIRED_REFIT_ROOT/SUCCESS"
(cd "$PAIRED_REFIT_ROOT"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Separate diagnostic refit bundle: $PAIRED_REFIT_ROOT"
