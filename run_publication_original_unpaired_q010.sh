#!/usr/bin/env bash
set -euo pipefail

# Authoritative end-to-end analysis for the manuscript's original, unpaired
# q <= 0.10 benchmark. This intentionally excludes the paired, disjoint, and
# feature-refitted exploratory workflows.

PROJECT="${PROJECT:-$PWD}"
RUN_ROOT="${RUN_ROOT:-RUNS_publication_original_unpaired_q010}"
R_BIN="${R_BIN:-}"
ALLOW_EXISTING="${ALLOW_EXISTING:-false}"
EXPECTED_R_VERSION="${EXPECTED_R_VERSION:-4.3.3}"
EXPECTED_MAASLIN2_VERSION="${EXPECTED_MAASLIN2_VERSION:-1.18.0}"

cd "$PROJECT"
PROJECT="$(pwd -P)"

if [[ -z "$R_BIN" ]]; then
  R_BIN="$(command -v Rscript || true)"
fi
if [[ -z "$R_BIN" || ! -x "$R_BIN" ]]; then
  echo "[ERROR] Rscript is unavailable. Set R_BIN to the container wrapper or an Rscript executable." >&2
  exit 1
fi

if [[ -e "$RUN_ROOT" && "$ALLOW_EXISTING" != "true" ]]; then
  echo "[ERROR] RUN_ROOT already exists: $RUN_ROOT" >&2
  echo "        Use a new RUN_ROOT. Set ALLOW_EXISTING=true only for an intentional resume." >&2
  exit 1
fi
mkdir -p "$RUN_ROOT/logs"
LOG="$RUN_ROOT/logs/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

required=(
  samples_spiked_all_independent.tsv
  samples_spiked_all_community.tsv
  spike_panel.tsv
  spike_taxon_aliases.csv
  metadata_w_study.tsv
  metaphlan4_merged_unspecified.csv
  kraken2_bracken_merged_unspecified.csv
  profile_results
  R/common_utils.R
  R/io_utils.R
  R/maaslin_utils.R
  R/spike_design_utils.R
  R/spike_metrics_utils.R
  R/spike_plot_utils.R
  R/spike_plot_utils_independent_patch.R
  R/classifier_utils.R
  scripts/00_build_spike_design.R
  scripts/01_compute_spike_metrics.R
  scripts/02_run_spike_biomarker_benchmark.R
  scripts/joint_correlate_01_02.R
  scripts/1st_panel_plot_manuscript_tool_discordance.R
  scripts/plot_manuscript_fnuc_case_study_panel.R
  scripts/plot_manuscript_independent_spike_overview.R
  scripts/plot_manuscript_crc_biomarker_recoverability_compact.R
  scripts/plot_manuscript_community_spike_concordance.R
  scripts/plot_manuscript_community_DA_recovery.R
  scripts/plot_internal_and_crossstudy_artefact_filter_validation.R
  scripts/plot_spike_calibratability_gam.R
  scripts/plot_manuscript_community_offtarget_artifacts.R
  scripts/plot_supp_community_depth_qc_recovery_models.R
  scripts/plot_supp_community_qc_recovery_full.R
  scripts/validate_taxon_aliases.R
  scripts/99_make_current_manuscript_figures.R
  rerun_supplementary_figures_only.sh
)
for path in "${required[@]}"; do
  [[ -e "$path" ]] || { echo "[ERROR] Missing required path: $path" >&2; exit 1; }
done

if grep -Eq 'make_option\("--random.effect"|random_effects[[:space:]]*=' scripts/02_run_spike_biomarker_benchmark.R; then
  echo "[ERROR] The primary 02 script appears to implement a random-effect/paired model." >&2
  exit 1
fi

echo "[INFO] Project:  $PROJECT"
echo "[INFO] Run root: $RUN_ROOT"
echo "[INFO] R runner: $R_BIN"
echo "[INFO] Started:  $(date --iso-8601=seconds)"

echo "[CHECK] Runtime and required R packages"
"$R_BIN" - "$EXPECTED_R_VERSION" "$EXPECTED_MAASLIN2_VERSION" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
expected_r <- args[[1]]
expected_maaslin <- args[[2]]
pkgs <- c("Maaslin2", "optparse", "data.table", "dplyr", "readr", "tidyr",
          "stringr", "tibble", "ggplot2", "forcats", "scales", "patchwork",
          "cowplot", "mgcv", "png")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing R packages: ", paste(pkgs[!ok], collapse = ", "))
actual_r <- paste(R.version$major, R.version$minor, sep = ".")
actual_maaslin <- as.character(utils::packageVersion("Maaslin2"))
if (actual_r != expected_r) stop("R version mismatch: expected ", expected_r, ", found ", actual_r)
if (actual_maaslin != expected_maaslin) {
  stop("MaAsLin2 version mismatch: expected ", expected_maaslin, ", found ", actual_maaslin)
}
cat("R: ", actual_r, "\nMaAsLin2: ", actual_maaslin, "\n", sep = "")
RS

echo "[CHECK] Complete taxon-alias mapping"
"$R_BIN" scripts/validate_taxon_aliases.R \
  --spike-panel spike_panel.tsv \
  --aliases spike_taxon_aliases.csv \
  --kraken kraken2_bracken_merged_unspecified.csv \
  --metaphlan metaphlan4_merged_unspecified.csv

echo "[STEP 1/8] Build spike design"
mkdir -p "$RUN_ROOT/design_auto"
"$R_BIN" scripts/00_build_spike_design.R \
  --independent_manifest samples_spiked_all_independent.tsv \
  --community_manifest samples_spiked_all_community.tsv \
  --spike_panel spike_panel.tsv \
  --outdir "$RUN_ROOT/design_auto"

echo "[STEP 2/8] Compute abundance-recovery metrics"
mkdir -p "$RUN_ROOT/spike_metrics"
"$R_BIN" scripts/01_compute_spike_metrics.R \
  --design "$RUN_ROOT/design_auto/spike_design.tsv" \
  --meta "$RUN_ROOT/design_auto/metadata_spiked.tsv" \
  --profile_root profile_results \
  --original_tables metaphlan4=metaphlan4_merged_unspecified.csv,kraken2_bracken=kraken2_bracken_merged_unspecified.csv \
  --auto_manifest_out "$RUN_ROOT/design_auto/run_manifest.auto.tsv" \
  --taxon_aliases spike_taxon_aliases.csv \
  --outdir "$RUN_ROOT/spike_metrics"

echo "[STEP 3/8] Run original unpaired MaAsLin2 benchmark"
mkdir -p "$RUN_ROOT/maaslin_spike"
"$R_BIN" scripts/02_run_spike_biomarker_benchmark.R \
  --design "$RUN_ROOT/design_auto/spike_design.tsv" \
  --meta "$RUN_ROOT/design_auto/metadata_spiked.tsv" \
  --metrics_dir "$RUN_ROOT/spike_metrics" \
  --run_manifest "$RUN_ROOT/design_auto/run_manifest.auto.tsv" \
  --filter_modes original \
  --background_conditions ALL \
  --fixed_effect spike_status \
  --reference_level unspiked \
  --q_threshold 0.10 \
  --outdir "$RUN_ROOT/maaslin_spike" \
  --taxon_aliases spike_taxon_aliases.csv

echo "[STEP 4/8] Compute biomarker-recoverability drivers"
mkdir -p "$RUN_ROOT/species_driver_and_thresholds"
"$R_BIN" scripts/joint_correlate_01_02.R \
  --metrics_dir "$RUN_ROOT/spike_metrics" \
  --maaslin_dir "$RUN_ROOT/maaslin_spike" \
  --outdir "$RUN_ROOT/species_driver_and_thresholds" \
  --focus_mode independent \
  --plot_filter_modes original \
  --detection_field positive \
  --min_detection_rate 0.5

echo "[STEP 5/8] Regenerate current manuscript Figures 2--5"
"$R_BIN" scripts/99_make_current_manuscript_figures.R \
  --project-root "$PROJECT" \
  --run-dir "$RUN_ROOT" \
  --metadata metadata_w_study.tsv \
  --kraken kraken2_bracken_merged_unspecified.csv \
  --metaphlan metaphlan4_merged_unspecified.csv

echo "[STEP 6/8] Generate paper Figure 6: post-hoc artefact exclusion"
mkdir -p "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion"
"$R_BIN" scripts/plot_internal_and_crossstudy_artefact_filter_validation.R \
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

echo "[STEP 7/8] Generate paper Figure 7: grouped-CV calibratability"
mkdir -p "$RUN_ROOT/manuscript_figures/figure7_calibratability"
"$R_BIN" scripts/plot_spike_calibratability_gam.R \
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

if [[ -n "${METRAPREP_ROOT:-}" || -n "${QC_LIST:-}" ||
      -n "${COMMUNITY_TARGET_FILE:-}" ||
      -s "$RUN_ROOT/qc_depth_original_samples/qc_files.txt" ]]; then
  echo "[SUPPLEMENT] Generate Supplementary Figures B1--B13"
  PROJECT="$PROJECT" \
  RUN_ROOT="$RUN_ROOT" \
  R_BIN="$R_BIN" \
  METRAPREP_ROOT="${METRAPREP_ROOT:-}" \
  QC_LIST="${QC_LIST:-}" \
  COMMUNITY_TARGET_FILE="${COMMUNITY_TARGET_FILE:-}" \
  bash rerun_supplementary_figures_only.sh
else
  echo "[WARN] Supplementary Figures B1--B13 were not generated because neither" >&2
  echo "       QC_LIST nor COMMUNITY_TARGET_FILE was supplied." >&2
  echo "       Run rerun_supplementary_figures_only.sh after supplying one of them." >&2
fi

echo "[STEP 8/8] Validate outputs and record provenance"
expected_outputs=(
  "$RUN_ROOT/manuscript_figures/current/Fig2_baseline_profiler_discordance.pdf"
  "$RUN_ROOT/manuscript_figures/current/Supplementary_Fig_B2_baseline_all_10_targets.pdf"
  "$RUN_ROOT/manuscript_figures/current/Fig3_independent_spike_recovery.pdf"
  "$RUN_ROOT/manuscript_figures/current/Fig4_biomarker_recoverability.pdf"
  "$RUN_ROOT/manuscript_figures/current/Fig5_community_spike_recovery_and_DA.pdf"
  "$RUN_ROOT/manuscript_figures/figure6_artefact_exclusion/manuscript_internal_and_crossstudy_artefact_filter_validation.pdf"
  "$RUN_ROOT/manuscript_figures/figure7_calibratability/manuscript_spike_calibratability_gam.pdf"
  "$RUN_ROOT/manuscript_figures/figure7_calibratability/grouped_cv_predictions.csv"
)
for path in "${expected_outputs[@]}"; do
  [[ -s "$path" ]] || { echo "[ERROR] Expected output is missing or empty: $path" >&2; exit 1; }
done

{
  echo "completed_at=$(date --iso-8601=seconds)"
  echo "project=$PROJECT"
  echo "run_root=$RUN_ROOT"
  echo "analysis=original_unpaired"
  echo "fixed_effect=spike_status"
  echo "reference=unspiked"
  echo "random_effect=none"
  echo "q_threshold=0.10"
  echo "calibratability_seed=1"
  "$R_BIN" -e 'cat("r_version=", paste(R.version$major,R.version$minor,sep="."), "\nMaaslin2_version=", as.character(packageVersion("Maaslin2")), "\n", sep="")'
  if command -v apptainer >/dev/null 2>&1 && [[ -n "${SIF:-}" && -s "${SIF:-}" ]]; then
    apptainer inspect "$SIF" 2>/dev/null || true
    sha256sum "$SIF"
  fi
} > "$RUN_ROOT/PROVENANCE.txt"

{
  for path in "${required[@]}"; do
    if [[ -f "$path" ]]; then
      sha256sum "$path"
    else
      find "$path" -type f -print0 | sort -z | xargs -0 sha256sum
    fi
  done
} > "$RUN_ROOT/INPUT_AND_CODE_SHA256.txt"
printf '%s\n' "${expected_outputs[@]}" > "$RUN_ROOT/EXPECTED_OUTPUTS.txt"

echo "[PASS] Complete original-unpaired paper analysis finished."
echo "[INFO] Output manifest: $RUN_ROOT/EXPECTED_OUTPUTS.txt"
echo "[INFO] Provenance:      $RUN_ROOT/PROVENANCE.txt"
echo "[INFO] Finished:        $(date --iso-8601=seconds)"
