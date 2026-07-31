#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
: "${SIF:?Set SIF to the MaAsLin2 analysis container image}"
R_BIN="${R_BIN:-}"
EXPECTED_R_VERSION="${EXPECTED_R_VERSION:-4.3.3}"
EXPECTED_MAASLIN2_VERSION="${EXPECTED_MAASLIN2_VERSION:-1.18.0}"

cd "$PROJECT"
fail=0
check_file() {
  if [[ ! -s "$1" ]]; then
    echo "[ERROR] Missing or empty file: $1" >&2
    fail=1
  else
    echo "[OK] $1"
  fi
}
check_dir() {
  if [[ ! -d "$1" ]] || ! find -L "$1" -type f -print -quit | grep -q .; then
    echo "[ERROR] Missing or empty directory: $1" >&2
    fail=1
  else
    echo "[OK] $1"
  fi
}

for file in \
  samples_spiked_all_independent.tsv samples_spiked_all_community.tsv \
  spike_panel.tsv spike_taxon_aliases.csv metadata_w_study.tsv \
  metaphlan4_merged_unspecified.csv kraken2_bracken_merged_unspecified.csv \
  run_publication_original_unpaired_q010.sh \
  run_original_unpaired_q010_cluster.sbatch \
  README.md \
  scripts/00_build_spike_design.R scripts/01_compute_spike_metrics.R \
  scripts/02_run_spike_biomarker_benchmark.R scripts/joint_correlate_01_02.R \
  scripts/1st_panel_plot_manuscript_tool_discordance.R \
  scripts/plot_manuscript_fnuc_case_study_panel.R \
  scripts/plot_manuscript_independent_spike_overview.R \
  scripts/plot_manuscript_crc_biomarker_recoverability_compact.R \
  scripts/plot_manuscript_community_spike_concordance.R \
  scripts/plot_manuscript_community_DA_recovery.R \
  scripts/plot_internal_and_crossstudy_artefact_filter_validation.R \
  scripts/plot_spike_calibratability_gam.R \
  scripts/plot_manuscript_community_offtarget_artifacts.R \
  scripts/plot_supplementary_taxon_specific_offtarget_patterns.R \
  scripts/plot_supp_community_depth_qc_recovery_models.R \
  scripts/plot_supp_community_qc_recovery_full.R \
  scripts/build_reference_representation_table.py \
  scripts/validate_taxon_aliases.R \
  scripts/99_make_current_manuscript_figures.R \
  prepare_reference_audit_assets.sh \
  mash_in_container.sh \
  rerun_current_figures_only.sh \
  rerun_supplementary_figures_only.sh
do
  check_file "$file"
done
for file in \
  R/common_utils.R R/io_utils.R R/maaslin_utils.R \
  R/spike_design_utils.R R/spike_metrics_utils.R R/spike_plot_utils.R \
  R/spike_plot_utils_independent_patch.R R/classifier_utils.R
do
  check_file "$file"
done
check_dir profile_results

if [[ "$fail" -ne 0 ]]; then
  echo "[FAIL] Required paths are incomplete; runtime checks were not attempted." >&2
  exit 1
fi

bash -n run_publication_original_unpaired_q010.sh \
  run_original_unpaired_q010_cluster.sbatch \
  cluster_preflight_original_unpaired.sh

if grep -Eq 'make_option\("--random.effect"|random_effects[[:space:]]*=' scripts/02_run_spike_biomarker_benchmark.R; then
  echo "[ERROR] Primary 02 script appears paired/random-effect based." >&2
  exit 1
fi

if [[ -z "$R_BIN" ]]; then
  check_file "$SIF"
  command -v apptainer >/dev/null 2>&1 || {
    echo "[ERROR] apptainer is unavailable." >&2
    exit 1
  }
  R_BIN="$PROJECT/rscript_in_container.sh"
fi
[[ -x "$R_BIN" ]] || { echo "[ERROR] R_BIN is not executable: $R_BIN" >&2; exit 1; }

PROJECT="$PROJECT" SIF="$SIF" "$R_BIN" - "$EXPECTED_R_VERSION" "$EXPECTED_MAASLIN2_VERSION" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
pkgs <- c("Maaslin2", "optparse", "data.table", "dplyr", "readr", "tidyr",
          "stringr", "tibble", "ggplot2", "forcats", "scales", "patchwork",
          "cowplot", "mgcv", "png")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Missing packages: ", paste(pkgs[!ok], collapse = ", "))
rv <- paste(R.version$major, R.version$minor, sep = ".")
mv <- as.character(packageVersion("Maaslin2"))
if (rv != args[[1]]) stop("Expected R ", args[[1]], "; found ", rv)
if (mv != args[[2]]) stop("Expected MaAsLin2 ", args[[2]], "; found ", mv)
cat("[OK] R ", rv, "; MaAsLin2 ", mv, "\n", sep = "")
RS

PROJECT="$PROJECT" SIF="$SIF" "$PROJECT/mash_in_container.sh" --version

PROJECT="$PROJECT" SIF="$SIF" "$R_BIN" - <<'RS'
files <- c(
  "samples_spiked_all_independent.tsv", "samples_spiked_all_community.tsv",
  "spike_panel.tsv", "spike_taxon_aliases.csv", "metadata_w_study.tsv",
  "metaphlan4_merged_unspecified.csv", "kraken2_bracken_merged_unspecified.csv"
)
for (f in files) {
  sep <- if (grepl("\\.tsv$", f)) "\t" else ","
  x <- utils::read.table(f, sep = sep, header = TRUE, nrows = 3,
                         quote = "\"", comment.char = "", check.names = FALSE)
  if (!ncol(x)) stop("No columns found in ", f)
  cat("[OK] readable: ", f, " (", ncol(x), " columns)\n", sep = "")
}
RS

PROJECT="$PROJECT" SIF="$SIF" "$R_BIN" - <<'RS'
helpers <- c(
  "R/common_utils.R", "R/io_utils.R", "R/maaslin_utils.R",
  "R/spike_design_utils.R", "R/spike_metrics_utils.R",
  "R/spike_plot_utils.R", "R/spike_plot_utils_independent_patch.R",
  "R/classifier_utils.R"
)
for (f in helpers) {
  tryCatch(
    parse(file = f),
    error = function(e) stop("Cannot parse required R helper ", f, ": ", conditionMessage(e))
  )
  cat("[OK] parseable: ", f, "\n", sep = "")
}
RS

echo
echo "[INFO] Validating complete 10-target x 2-profiler alias mapping..."
PROJECT="$PROJECT" SIF="$SIF" "$R_BIN" scripts/validate_taxon_aliases.R \
  --spike-panel spike_panel.tsv \
  --aliases spike_taxon_aliases.csv \
  --kraken kraken2_bracken_merged_unspecified.csv \
  --metaphlan metaphlan4_merged_unspecified.csv

echo "[PASS] Preflight passed."
