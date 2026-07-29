#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
: "${RUN_ROOT:?Set RUN_ROOT to the completed analysis directory}"
R_BIN="${R_BIN:-$PROJECT/rscript_in_container.sh}"
QC_LIST="${QC_LIST:-}"
COMMUNITY_TARGET_FILE="${COMMUNITY_TARGET_FILE:-}"
METAPREP_ROOT="${METAPREP_ROOT:-}"
METADATA="${METADATA:-}"

cd "$PROJECT"
PROJECT="$(pwd -P)"
[[ "$RUN_ROOT" = /* ]] || RUN_ROOT="$PROJECT/$RUN_ROOT"
METADATA="${METADATA:-$PROJECT/metadata_w_study.tsv}"
[[ "$METADATA" = /* ]] || METADATA="$PROJECT/$METADATA"

[[ -d "$RUN_ROOT" ]] || {
  echo "[ERROR] Completed RUN_ROOT does not exist: $RUN_ROOT" >&2
  exit 1
}
[[ -x "$R_BIN" ]] || {
  echo "[ERROR] R_BIN is not executable: $R_BIN" >&2
  exit 1
}

current="$RUN_ROOT/manuscript_figures/current"
panels="$current/source_panels"
outdir="$RUN_ROOT/manuscript_figures/supplementary"
source="$outdir/source_panels"
qc_source="$source/community_qc_recovery"
offtarget_source="$source/community_offtarget_artifacts"
baseline_source="$source/baseline_all_10_targets"
independent_source="$source/independent_all_fractions"
community_da_source="$source/community_DA_all_fractions"
mkdir -p \
  "$outdir" "$source" "$qc_source" "$offtarget_source" \
  "$baseline_source" "$independent_source" "$community_da_source"

run_r() {
  "$R_BIN" "$@"
}

copy_pair() {
  local input_stem="$1"
  local output_stem="$2"
  local ext
  for ext in pdf png; do
    [[ -s "$input_stem.$ext" ]] || {
      echo "[ERROR] Missing supplementary source: $input_stem.$ext" >&2
      exit 1
    }
    cp -f "$input_stem.$ext" "$outdir/$output_stem.$ext"
  done
  echo "[OK] $output_stem"
}

echo "[INFO] Reusing completed analyses in: $RUN_ROOT"
echo "[INFO] No MaAsLin2 or abundance-recovery analysis will be rerun."

# Regenerate, rather than copy, all panels whose plot scripts may have changed.
# This keeps a supplementary-only rerun faithful to the current repository.
run_r "$PROJECT/scripts/1st_panel_plot_manuscript_tool_discordance.R" \
  --indir "$RUN_ROOT/spike_metrics" \
  --outdir "$baseline_source" \
  --spike-labels "Bfrag,Csym,Dpne,Fnuc,Hhat,Pmic,Pana,Psto,Porp,Pint" \
  --metadata-input "$METADATA" \
  --metadata-sample-col sample_id \
  --metadata-condition-col Target_Condition \
  --metadata-study-col Study \
  --width 12.5 \
  --height 16.5 \
  --dpi 450

run_r "$PROJECT/scripts/plot_manuscript_independent_spike_overview.R" \
  --indir "$RUN_ROOT/spike_metrics" \
  --outdir "$independent_source" \
  --spike-labels "Bfrag,Csym,Dpne,Fnuc,Hhat,Pmic,Pana,Psto,Porp,Pint" \
  --sort-labels-by taxon \
  --dpi 450

run_r "$PROJECT/scripts/plot_manuscript_community_DA_recovery.R" \
  --indir "$RUN_ROOT/maaslin_spike" \
  --outdir "$community_da_source" \
  --community-size 10 \
  --main-effective-fractions "0.00001,0.00005,0.0001,0.0005,0.001" \
  --spike-label-order "Bfrag,Csym,Dpne,Fnuc,Hhat,Pmic,Pana,Psto,Porp,Pint" \
  --dpi 450

# B1 and B6--B10 require original-sample read depth. Reuse the derived table
# when available; otherwise build it from the FastQC zip-file inventory.
if [[ -z "$COMMUNITY_TARGET_FILE" ]]; then
  candidates=(
    "$RUN_ROOT/manuscript_figures/supplementary/source_panels/community_depth_qc/community_target_level_depth_recovery.tsv"
    "$RUN_ROOT/revised_manuscript_figures/source_panels/FigC8_community_depth_qc_recovery/community_target_level_depth_recovery.tsv"
  )
  for candidate in "${candidates[@]}"; do
    if [[ -s "$candidate" ]]; then
      COMMUNITY_TARGET_FILE="$candidate"
      break
    fi
  done
fi

if [[ -z "$COMMUNITY_TARGET_FILE" ]]; then
  if [[ -z "$QC_LIST" && -s "$RUN_ROOT/qc_depth_original_samples/qc_files.txt" ]]; then
    QC_LIST="$RUN_ROOT/qc_depth_original_samples/qc_files.txt"
  fi
  if [[ -z "$QC_LIST" && -n "$METAPREP_ROOT" ]]; then
    qc_inventory="$RUN_ROOT/qc_depth_original_samples"
    mkdir -p "$qc_inventory"
    python3 "$PROJECT/scripts/build_metaprep_qc_manifest.py" \
      --metaprep-root "$METAPREP_ROOT" \
      --metadata "$METADATA" \
      --output-list "$qc_inventory/qc_files.txt" \
      --output-manifest "$qc_inventory/qc_manifest.tsv"
    QC_LIST="$qc_inventory/qc_files.txt"
  fi
  [[ -n "$QC_LIST" && -s "$QC_LIST" ]] || {
    cat >&2 <<EOF
[ERROR] Supplementary Figures B1 and B6--B10 require original-sample QC depth.
Set one of:
  COMMUNITY_TARGET_FILE=/path/to/community_target_level_depth_recovery.tsv
or:
  METAPREP_ROOT=/path/to/data/processed/metaprep
or:
  QC_LIST=/path/to/qc_files.txt

qc_files.txt must contain one absolute path per FastQC *_fastqc.zip file from
the pre- and post-QC original samples.
EOF
    exit 1
  }
  qc_build="$source/community_depth_qc"
  mkdir -p "$qc_build"
  run_r "$PROJECT/scripts/plot_supp_community_depth_qc_recovery_models.R" \
    --run-root "$RUN_ROOT" \
    --metadata "$PROJECT/metadata_w_study.tsv" \
    --qc-list "$QC_LIST" \
    --recovery-file "$RUN_ROOT/spike_metrics/target_member_errors_with_condition.csv" \
    --outdir "$qc_build"
  COMMUNITY_TARGET_FILE="$qc_build/community_target_level_depth_recovery.tsv"
fi

[[ -s "$COMMUNITY_TARGET_FILE" ]] || {
  echo "[ERROR] Missing community/QC recovery table: $COMMUNITY_TARGET_FILE" >&2
  exit 1
}

run_r "$PROJECT/scripts/plot_supp_community_qc_recovery_full.R" \
  --run-root "$RUN_ROOT" \
  --community-target-file "$COMMUNITY_TARGET_FILE" \
  --outdir "$qc_source" \
  --dpi 450

# B13 is a plotting-only analysis of completed abundance-error and MaAsLin2
# outputs.
run_r "$PROJECT/scripts/plot_manuscript_community_offtarget_artifacts.R" \
  --indir "$RUN_ROOT/maaslin_spike" \
  --outdir "$offtarget_source" \
  --alias-file "$PROJECT/spike_taxon_aliases.csv" \
  --community-size 10 \
  --filter-mode original \
  --dpi 450

# Collect the freshly rendered, publication-sized panels under stable B1--B13
# names. PDFs are the publication masters; PNGs are high-resolution previews.
copy_pair \
  "$qc_source/FigB8_QC_only_ABC" \
  "FigB1_sequencing_depth_and_preprocessing_QC"
copy_pair \
  "$baseline_source/manuscript_full_baseline_discordance_panel" \
  "FigB2_baseline_all_10_targets"
copy_pair \
  "$independent_source/supp_good_recovery_dot_heatmap_allfractions" \
  "FigB3_independent_good_recovery_all_fractions"
copy_pair \
  "$independent_source/supp_recovery_class_composition_allfractions" \
  "FigB4_independent_recovery_classes_all_fractions"
copy_pair \
  "$independent_source/supp_bias_variability_dotplots_allfractions" \
  "FigB5_independent_bias_variability_all_fractions"
copy_pair \
  "$qc_source/FigB9_recovery_classes_vs_spike_fraction" \
  "FigB6_community_recovery_classes"
copy_pair \
  "$qc_source/FigB11_taxon_good_recovery_vs_spike_fraction" \
  "FigB7_community_taxon_good_recovery"
copy_pair \
  "$qc_source/FigB12_taxon_average_intermediate_recovery_vs_spike_fraction" \
  "FigB8_community_taxon_intermediate_recovery"
copy_pair \
  "$qc_source/FigB13_taxon_poor_missed_recovery_vs_spike_fraction" \
  "FigB9_community_taxon_poor_missed_recovery"
copy_pair \
  "$qc_source/FigB10_recovery_class_composition_by_read_support" \
  "FigB10_community_recovery_by_read_support"
copy_pair \
  "$community_da_source/supp_panel_B_community_target_DA_detection_all_effective_fractions" \
  "FigB11_community_target_DA_all_fractions"
copy_pair \
  "$community_da_source/supp_panel_D_community_offtarget_enriched_DA_burden_all_effective_fractions" \
  "FigB12_community_offtarget_DA_burden"
copy_pair \
  "$offtarget_source/manuscript_community_offtarget_artifacts_overview" \
  "FigB13_offtarget_DA_profiler_artefacts"

{
  echo "Supplementary figures regenerated from completed run: $RUN_ROOT"
  echo "Community/QC input: $COMMUNITY_TARGET_FILE"
  find "$outdir" -maxdepth 1 -type f \( -name 'FigB*.pdf' -o -name 'FigB*.png' \) \
    -printf '%f\n' | sort
} > "$outdir/SUPPLEMENTARY_FIGURE_MANIFEST.txt"

pdf_count="$(find "$outdir" -maxdepth 1 -type f -name 'FigB*.pdf' | wc -l)"
png_count="$(find "$outdir" -maxdepth 1 -type f -name 'FigB*.png' | wc -l)"
if [[ "$pdf_count" -ne 13 || "$png_count" -ne 13 ]]; then
  echo "[ERROR] Expected 13 PDFs and 13 PNGs; found $pdf_count PDFs and $png_count PNGs." >&2
  exit 1
fi

# Reject all-white or near-blank review images before reporting success. This
# catches the graphics-device failure that previously produced empty current
# figures even though their source panels were valid.
run_r - "$outdir" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[[1]]
files <- sort(list.files(outdir, pattern = "^FigB[0-9]+.*\\.png$", full.names = TRUE))
if (!requireNamespace("png", quietly = TRUE)) {
  stop("The R package 'png' is required for supplementary image validation.")
}
bad <- character()
for (f in files) {
  img <- png::readPNG(f)
  rgb <- img[, , seq_len(min(3L, dim(img)[3L])), drop = FALSE]
  # Inspect at most a 600 x 600 representative grid. Applying a function to
  # every pixel of a 450-dpi multi-panel image is needlessly slow and can look
  # like a hung process on a login node.
  row_idx <- unique(round(seq(1, dim(rgb)[1], length.out = min(600L, dim(rgb)[1]))))
  col_idx <- unique(round(seq(1, dim(rgb)[2], length.out = min(600L, dim(rgb)[2]))))
  sample_rgb <- rgb[row_idx, col_idx, , drop = FALSE]
  pixel_min <- do.call(pmin, lapply(seq_len(dim(sample_rgb)[3]), function(k) sample_rgb[, , k]))
  nonwhite <- mean(pixel_min < 0.985)
  spread <- diff(range(sample_rgb, finite = TRUE))
  if (!is.finite(nonwhite) || nonwhite < 0.002 || spread < 0.05) bad <- c(bad, basename(f))
}
if (length(bad)) stop("Blank or near-blank supplementary PNG(s): ", paste(bad, collapse = ", "))
message("[PASS] All 13 supplementary PNGs contain non-white graphical content.")
RS

echo "[PASS] Supplementary Figures B1--B13 regenerated."
echo "[INFO] Output: $outdir"
