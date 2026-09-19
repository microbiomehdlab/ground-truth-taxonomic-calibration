#!/usr/bin/env bash
set -euo pipefail

: "${PROFILE_MANIFEST:?Set PROFILE_MANIFEST}"
: "${ENDPOINTS:?Set ENDPOINTS}"
: "${ABUNDANCE:?Set ABUNDANCE}"
: "${OUTDIR:?Set OUTDIR}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STATUS="${REPORT_STATUS:-DEVELOPMENT_ONLY}"
mkdir -p "$OUTDIR"

builder_args=()
if [[ -n "${FEATURE_ALIASES:-}" ]]; then
  builder_args+=(--feature-aliases "$FEATURE_ALIASES")
fi
python3 "$ROOT/analysis_v2/scripts/build_perturbation_response_input.py" \
  --profile-manifest "$PROFILE_MANIFEST" \
  --endpoints "$ENDPOINTS" \
  --abundance "$ABUNDANCE" \
  "${builder_args[@]}" \
  --outdir "$OUTDIR/input"

analysis_args=(
  --input "$OUTDIR/input/paired_feature_responses.parquet"
  --outdir "$OUTDIR/analysis"
)
if [[ -n "${DISEASE_RESULTS:-}" ]]; then
  analysis_args+=(--disease-results "$DISEASE_RESULTS")
fi
if [[ -n "${FEATURE_ALIASES:-}" ]]; then
  analysis_args+=(--feature-aliases "$FEATURE_ALIASES")
fi
python3 "$ROOT/analysis_v2/scripts/analyze_perturbation_response.py" "${analysis_args[@]}"

recovery_args=(
  --responses "$OUTDIR/input/paired_feature_responses.parquet"
  --outdir "$OUTDIR/target_recovery"
)
if [[ -n "${BIOMARKER_METRICS:-}" ]]; then
  recovery_args+=(--biomarker-metrics "$BIOMARKER_METRICS")
fi
python3 "$ROOT/analysis_v2/scripts/summarize_target_recovery.py" "${recovery_args[@]}"

Rscript "$ROOT/analysis_v2/scripts/make_perturbation_response_report.R" \
  --input-root "$OUTDIR/analysis" \
  --target-recovery-root "$OUTDIR/target_recovery" \
  --outdir "$OUTDIR/report" \
  --report-status "$STATUS"

for required in \
  input/SUCCESS analysis/SUCCESS analysis/perturbation_response.sha256 \
  target_recovery/SUCCESS target_recovery/target_recovery.sha256 \
  analysis/response_operator.tsv analysis/reliability_certificates.tsv \
  analysis/heldout_validation.tsv analysis/panel_saturation.tsv \
  report/SUCCESS report/provenance/perturbation_response_report.sha256
do
  test -s "$OUTDIR/$required" || { echo "[ERROR] Missing $required" >&2; exit 1; }
done

printf 'status=%s\nuse_for_manuscript=%s\n' "$STATUS" \
  "$([[ "$STATUS" == "DEVELOPMENT_ONLY" ]] && echo NO || echo YES)" \
  > "$OUTDIR/DEVELOPMENT_ONLY.txt"
find "$OUTDIR/input" "$OUTDIR/analysis" "$OUTDIR/target_recovery" "$OUTDIR/report" \
  -type f -print0 |
  sort -z | xargs -0 sha256sum > "$OUTDIR/perturbation_response_atlas.sha256"
printf 'status\tPASS\nreport_status\t%s\n' "$STATUS" > "$OUTDIR/SUCCESS"
echo "[PASS] Perturbation-response atlas completed: $OUTDIR"
