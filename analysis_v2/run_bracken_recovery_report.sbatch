#!/bin/bash
#SBATCH --job-name=bracken_report
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=bracken_report.%j.out
#SBATCH --error=bracken_report.%j.err
set -euo pipefail
: "${PROJECT:?Export repository root}"
: "${RECOVERY_ROOT:?Export the verified three-cohort recovery root}"
: "${REPORT_ROOT:?Export a fresh output root}"
: "${ANALYSIS_SIF:?Export the pinned analysis image}"
cd "$PROJECT"
test ! -e "$REPORT_ROOT" || { echo 'Choose a fresh REPORT_ROOT' >&2; exit 1; }
test -s "$ANALYSIS_SIF"
appt=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF")
"${appt[@]}" Rscript -e 'stopifnot(requireNamespace("ggplot2",quietly=TRUE),requireNamespace("scales",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE))'
python3 -B analysis_v2/scripts/summarize_bracken_recovery.py \
    --recovery-root "$RECOVERY_ROOT" --outdir "$REPORT_ROOT/tables"
"${appt[@]}" Rscript analysis_v2/scripts/plot_bracken_recovery_comparison.R \
    --summary-dir "$REPORT_ROOT/tables" --outdir "$REPORT_ROOT/figures"
git rev-parse HEAD > "$REPORT_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$REPORT_ROOT/analysis_image.sha256"
sha256sum analysis_v2/scripts/summarize_bracken_recovery.py \
    analysis_v2/scripts/plot_bracken_recovery_comparison.R > "$REPORT_ROOT/report_code.sha256"
(
    cd "$REPORT_ROOT"
    find . -type f ! -name SHA256SUMS ! -name SUCCESS -print0 |
        sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c --quiet SHA256SUMS
)
printf 'status\tPASS_REPORT\n' > "$REPORT_ROOT/SUCCESS"
echo "[PASS] Three-cohort Bracken report: $REPORT_ROOT"
