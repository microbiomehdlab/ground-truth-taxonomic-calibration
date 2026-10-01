#!/bin/bash
#SBATCH --job-name=metaphlan_report
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=metaphlan_report.%j.out
#SBATCH --error=metaphlan_report.%j.err
set -euo pipefail
: "${PROJECT:?Export project root}"
: "${METAPHLAN_RECOVERY_ROOT:?Export verified endpoint bundle}"
: "${MPA_REPORT_ROOT:?Export fresh report root}"
: "${ANALYSIS_SIF:?Export pinned image}"
cd "$PROJECT"
test ! -e "$MPA_REPORT_ROOT"
appt=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF")
"${appt[@]}" Rscript -e 'stopifnot(requireNamespace("ggplot2",quietly=TRUE),requireNamespace("scales",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE))'
"${appt[@]}" python3 -B analysis_v2/scripts/summarize_metaphlan_recovery.py \
  --recovery-root "$METAPHLAN_RECOVERY_ROOT" --outdir "$MPA_REPORT_ROOT/tables"
"${appt[@]}" Rscript analysis_v2/scripts/plot_bracken_recovery_comparison.R \
  --summary-dir "$MPA_REPORT_ROOT/tables" --outdir "$MPA_REPORT_ROOT/figures"
git rev-parse HEAD > "$MPA_REPORT_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$MPA_REPORT_ROOT/analysis_image.sha256"
sha256sum analysis_v2/scripts/summarize_metaphlan_recovery.py \
  analysis_v2/scripts/summarize_bracken_recovery.py \
  analysis_v2/scripts/plot_bracken_recovery_comparison.R > "$MPA_REPORT_ROOT/report_code.sha256"
printf 'status\tPASS_REPORT\nprofiler\tmetaphlan4\n' > "$MPA_REPORT_ROOT/SUCCESS"
(cd "$MPA_REPORT_ROOT"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] MetaPhlAn reference-comparison report: $MPA_REPORT_ROOT"
