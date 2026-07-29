#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
: "${RUN_ROOT:?Set RUN_ROOT to the completed analysis directory}"
R_BIN="${R_BIN:-$PROJECT/rscript_in_container.sh}"
RERUN_PANELS="${RERUN_PANELS:-false}"
RERUN_FIGURES_6_7="${RERUN_FIGURES_6_7:-true}"

cd "$PROJECT"
PROJECT="$(pwd -P)"

if [[ "$RUN_ROOT" != /* ]]; then
  RUN_ROOT="$PROJECT/$RUN_ROOT"
fi
[[ -d "$RUN_ROOT" ]] || {
  echo "[ERROR] Completed RUN_ROOT does not exist: $RUN_ROOT" >&2
  exit 1
}
[[ -x "$R_BIN" ]] || {
  echo "[ERROR] R_BIN is not executable: $R_BIN" >&2
  exit 1
}

required=(
  "$PROJECT/scripts/99_make_current_manuscript_figures.R"
)
if [[ "$RERUN_PANELS" != "true" ]]; then
  required+=("$RUN_ROOT/manuscript_figures/current/source_panels")
fi
if [[ "$RERUN_FIGURES_6_7" == "true" ]]; then
  required+=(
    "$PROJECT/scripts/plot_internal_and_crossstudy_artefact_filter_validation.R"
    "$PROJECT/scripts/plot_spike_calibratability_gam.R"
    "$RUN_ROOT/spike_metrics/community/species_trace_with_condition.csv"
    "$RUN_ROOT/maaslin_spike/maaslin_significant_features_ALLFILTERS.csv"
    "$RUN_ROOT/spike_metrics/community/target_member_errors_with_condition.csv"
  )
fi
for path in "${required[@]}"; do
  [[ -e "$path" ]] || {
    echo "[ERROR] Missing required path: $path" >&2
    exit 1
  }
done

backup="$RUN_ROOT/manuscript_figures/previous_main_figures_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$backup"
for stem in \
  Fig2_baseline_profiler_discordance \
  Supplementary_Fig_B2_baseline_all_10_targets \
  Fig3_independent_spike_recovery \
  Fig4_biomarker_recoverability \
  Fig5_community_spike_recovery_and_DA
do
  for ext in pdf png; do
    path="$RUN_ROOT/manuscript_figures/current/$stem.$ext"
    [[ -f "$path" ]] && mv "$path" "$backup/"
  done
done
if [[ "$RERUN_FIGURES_6_7" == "true" ]]; then
  for path in \
    "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion/manuscript_internal_and_crossstudy_artefact_filter_validation.pdf" \
    "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion/manuscript_internal_and_crossstudy_artefact_filter_validation.png" \
    "$RUN_ROOT/manuscript_figures/figure7_calibratability/manuscript_spike_calibratability_gam.pdf" \
    "$RUN_ROOT/manuscript_figures/figure7_calibratability/manuscript_spike_calibratability_gam.png"
  do
    [[ -f "$path" ]] && mv "$path" "$backup/"
  done
fi

echo "[INFO] Project:             $PROJECT"
echo "[INFO] Completed run:       $RUN_ROOT"
echo "[INFO] Previous figures:    $backup"
if [[ "$RERUN_PANELS" == "true" ]]; then
  echo "[INFO] Regenerating source panels from completed analysis tables; statistical analyses will not be rerun."
else
  echo "[INFO] Reusing source panels; analyses will not be rerun."
fi
if [[ "$RERUN_FIGURES_6_7" == "true" ]]; then
  echo "[INFO] Figures 6--7 will be regenerated from completed upstream tables."
  echo "[INFO] Figure 7 repeats its deterministic, seeded grouped-CV plotting fit."
fi

"$R_BIN" "$PROJECT/scripts/99_make_current_manuscript_figures.R" \
  --project-root "$PROJECT" \
  --run-dir "$RUN_ROOT" \
  --metadata metadata_w_study.tsv \
  --kraken kraken2_bracken_merged_unspecified.csv \
  --metaphlan metaphlan4_merged_unspecified.csv \
  --rerun-panels "$RERUN_PANELS"

if [[ "$RERUN_FIGURES_6_7" == "true" ]]; then
  mkdir -p "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion"
  "$R_BIN" "$PROJECT/scripts/plot_internal_and_crossstudy_artefact_filter_validation.R" \
    --trace "$RUN_ROOT/spike_metrics/community/species_trace_with_condition.csv" \
    --significant-file "$RUN_ROOT/maaslin_spike/maaslin_significant_features_ALLFILTERS.csv" \
    --outdir "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion" \
    --q-threshold 0.10 \
    --threshold 0.05 \
    --thresholds 0.025,0.05,0.10,0.20 \
    --min-samples 20 \
    --min-median-expected 0 \
    --community-label CRCpanel \
    --filter-mode original

  mkdir -p "$RUN_ROOT/manuscript_figures/figure7_calibratability"
  "$R_BIN" "$PROJECT/scripts/plot_spike_calibratability_gam.R" \
    --input "$RUN_ROOT/spike_metrics/community/target_member_errors_with_condition.csv" \
    --outdir "$RUN_ROOT/manuscript_figures/figure7_calibratability" \
    --folds 5 \
    --repeats 5 \
    --seed 1 \
    --k 5 \
    --min-rows 80 \
    --min-detected 40 \
    --focus-targets "Bacteroides fragilis,Fusobacterium nucleatum subsp. nucleatum,Parvimonas micra,Dialister pneumosintes" \
    --transfer-tests true

  "$R_BIN" - \
    "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion/manuscript_internal_and_crossstudy_artefact_filter_validation.png" \
    "$RUN_ROOT/manuscript_figures/figure7_calibratability/manuscript_spike_calibratability_gam.png" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
if (!requireNamespace("png", quietly = TRUE)) {
  stop("The R package 'png' is required for final-figure validation.")
}
for (path in args) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing or empty final figure: ", path)
  }
  image <- png::readPNG(path)
  rgb <- image[, , seq_len(min(3L, dim(image)[3L])), drop = FALSE]
  row_idx <- unique(round(seq(1, dim(rgb)[1], length.out = min(600L, dim(rgb)[1]))))
  col_idx <- unique(round(seq(1, dim(rgb)[2], length.out = min(600L, dim(rgb)[2]))))
  sample_rgb <- rgb[row_idx, col_idx, , drop = FALSE]
  nonwhite <- mean(sample_rgb < 0.985, na.rm = TRUE)
  spread <- diff(range(sample_rgb, finite = TRUE))
  if (!is.finite(nonwhite) || nonwhite < 0.001 ||
      !is.finite(spread) || spread < 0.02) {
    stop("Blank or near-blank final figure: ", path)
  }
  message("[OK] Non-blank final figure: ", basename(path))
}
RS
fi

echo "[PASS] Manuscript Figures 2--7 regenerated from completed analysis outputs."
echo "[INFO] Figures 2--5: $RUN_ROOT/manuscript_figures/current"
if [[ "$RERUN_FIGURES_6_7" == "true" ]]; then
  echo "[INFO] Figure 6: $RUN_ROOT/manuscript_figures/figure6_artefact_exclusion"
  echo "[INFO] Figure 7: $RUN_ROOT/manuscript_figures/figure7_calibratability"
fi
