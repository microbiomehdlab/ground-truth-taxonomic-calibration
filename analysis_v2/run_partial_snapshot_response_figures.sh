#!/usr/bin/env bash
# Development-only response figures from a frozen, partially completed
# three-cohort snapshot. Never uses the older read-reference defaults.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
: "${CANONICAL_INPUT:?Set CANONICAL_INPUT to the frozen three-cohort canonical input}"
: "${ENDPOINTS:?Set ENDPOINTS to corrected three-cohort paired_endpoints.tsv}"
: "${TARGET_GENOME_SIZES:?Set TARGET_GENOME_SIZES to measured spike FASTA lengths}"
: "${EFFECTIVE_GENOME_SIZES:?Set EFFECTIVE_GENOME_SIZES to audited three-cohort G_eff}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF to the validated analysis image}"
: "${EXPECTED_BASELINES:?Set EXPECTED_BASELINES to the snapshot sample-wide count}"
: "${OUTDIR:?Set OUTDIR to a new development-only output directory}"
ALIASES="${ALIASES:-examples/spike_taxon_aliases.csv}"
DISEASE_CALLS=""
if [[ -n "${SOURCE_RUN:-}" ]]; then
  DISEASE_CALLS="$SOURCE_RUN/models/disease/models/primary_disease_da_results.tsv"
  for path in "$SOURCE_RUN/SUCCESS" "$SOURCE_RUN/DEVELOPMENT_ONLY.txt" \
              "$SOURCE_RUN/models/disease/SUCCESS" "$DISEASE_CALLS"; do
    [[ -s "$path" ]] || { echo "[ERROR] Missing disease source: $path" >&2; exit 1; }
  done
fi
for path in "$CANONICAL_INPUT" "$ENDPOINTS" "$TARGET_GENOME_SIZES" \
            "$EFFECTIVE_GENOME_SIZES" "$ANALYSIS_SIF" "$ALIASES"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing input: $path" >&2; exit 1; }
done
[[ "$EXPECTED_BASELINES" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] EXPECTED_BASELINES must be a positive integer" >&2; exit 1;
}
[[ ! -e "$OUTDIR" ]] || { echo "[ERROR] OUTDIR exists: $OUTDIR" >&2; exit 1; }
mkdir -p "$OUTDIR"
printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"

analysis_python() {
  apptainer exec --cleanenv --bind "$ROOT:$ROOT" --pwd "$ROOT" \
    "$ANALYSIS_SIF" python3 "$@"
}
analysis_r() {
  apptainer exec --cleanenv --bind "$ROOT:$ROOT" --pwd "$ROOT" \
    "$ANALYSIS_SIF" Rscript "$@"
}

# The baseline manifest is the frozen unit inventory for this provisional run.
analysis_python analysis_v2/scripts/select_baseline_manifest.py \
  --canonical "$CANONICAL_INPUT" --expected-profiles "$EXPECTED_BASELINES" \
  --outdir "$OUTDIR/baseline_selection"
test -s "$OUTDIR/baseline_selection/SUCCESS"
awk -F '\t' '
  NR == 1 { for (i=1; i<=NF; i++) if ($i == "cohort") col=i; next }
  { n[$col]++ }
  END {
    for (k in n) count++
    if (count != 3 || !n["feng"] || !n["yachida"] || !n["zeller"])
      exit 1
    printf "[INFO] snapshot baselines: feng=%d yachida=%d zeller=%d\n",
           n["feng"], n["yachida"], n["zeller"]
  }' "$OUTDIR/baseline_selection/baseline_manifest.tsv" || {
    echo "[ERROR] snapshot is not exactly Feng/Yachida/Zeller" >&2; exit 1;
  }
analysis_python analysis_v2/scripts/build_biomarker_abundance_input.py \
  --canonical "$CANONICAL_INPUT" --aliases "$ALIASES" \
  --outdir "$OUTDIR/native_abundance"
test -s "$OUTDIR/native_abundance/SUCCESS"

analysis_python analysis_v2/scripts/build_perturbation_response_input.py \
  --profile-manifest "$OUTDIR/native_abundance/biomarker_profile_manifest.tsv" \
  --endpoints "$ENDPOINTS" \
  --abundance "$OUTDIR/native_abundance/biomarker_abundance_long.tsv" \
  --metaphlan-reference genome_equivalent \
  --target-genome-sizes "$TARGET_GENOME_SIZES" \
  --effective-genome-sizes "$EFFECTIVE_GENOME_SIZES" \
  --outdir "$OUTDIR/response_input"
test -s "$OUTDIR/response_input/SUCCESS"
RESPONSES="$OUTDIR/response_input/paired_feature_responses.parquet"
analyzer_args=(--responses "$RESPONSES" --reference-scale profiler_scale
  --cohort-validation require_holdout --feature-aliases "$ALIASES"
  --outdir "$OUTDIR/response_analysis")
if [[ -n "$DISEASE_CALLS" ]]; then
  analyzer_args+=(--disease-results "$DISEASE_CALLS")
fi
analysis_python analysis_v2/scripts/analyze_perturbation_response.py "${analyzer_args[@]}"
analysis_python analysis_v2/scripts/summarize_target_recovery.py \
  --responses "$RESPONSES" --reference-scale profiler_scale \
  --outdir "$OUTDIR/target_recovery"
analysis_r analysis_v2/scripts/make_perturbation_response_report.R \
  --input-root "$OUTDIR/response_analysis" \
  --target-recovery-root "$OUTDIR/target_recovery" \
  --report-status DEVELOPMENT_ONLY \
  --outdir "$OUTDIR/report"
for path in "$OUTDIR/response_analysis/SUCCESS" \
            "$OUTDIR/target_recovery/SUCCESS" "$OUTDIR/report/SUCCESS" \
            "$OUTDIR/report/figures/response_operator_crosstalk.png" \
            "$OUTDIR/report/figures/response_operator_performance_map.png" \
            "$OUTDIR/report/figures/superposition_error_comparison.png" \
            "$OUTDIR/report/figures/superposition_improvement_by_dose.png"; do
  test -s "$path" || { echo "[ERROR] Missing output: $path" >&2; exit 1; }
done
if [[ -n "$DISEASE_CALLS" ]]; then
  test -s "$OUTDIR/report/figures/biomarker_replication_by_reliability.png" || {
    echo "[ERROR] Missing external-replication candidate plot" >&2; exit 1;
  }
fi
sha256sum "$CANONICAL_INPUT" "$ENDPOINTS" "$TARGET_GENOME_SIZES" \
  "$EFFECTIVE_GENOME_SIZES" "$ALIASES" "$ANALYSIS_SIF" \
  > "$OUTDIR/run_inputs.sha256"
if [[ -n "$DISEASE_CALLS" ]]; then
  sha256sum "$SOURCE_RUN/SUCCESS" "$DISEASE_CALLS" >> "$OUTDIR/run_inputs.sha256"
fi
printf 'status\tPASS\nanalysis_status\tDEVELOPMENT_ONLY\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Partial-snapshot profiler-scale response figures: $OUTDIR"
