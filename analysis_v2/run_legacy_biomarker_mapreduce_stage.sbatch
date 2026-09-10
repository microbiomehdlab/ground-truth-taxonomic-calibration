#!/usr/bin/env bash
set -euo pipefail
ROOT="${PROJECT:-${SLURM_SUBMIT_DIR:-}}"; cd "$ROOT"
: "${PIPELINE_STAGE:?}"; : "${DEV_INPUT_ROOT:?}"; : "${ANALYSIS_SIF:?}"; : "${OUTDIR:?}"
canonical="$DEV_INPUT_ROOT/canonical_input.tsv"; metadata="$OUTDIR/metadata/sample_metadata.tsv"
inputs="$OUTDIR/model_inputs"; aliases=examples/spike_taxon_aliases.csv; panel=spikes/spike_panel.tsv
case "$PIPELINE_STAGE" in
 setup)
  mkdir -p "$OUTDIR"/{metadata,models,reports,provenance}
  printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt"
  python3 analysis_v2/scripts/build_development_sample_metadata.py --canonical "$canonical" \
    --feng-manifest datasets/fengq/manifests/production_manifest.tsv --zeller-manifest datasets/zellerg/manifests/production_manifest.tsv --output "$metadata"
  python3 analysis_v2/scripts/build_biomarker_abundance_input.py --canonical "$canonical" --aliases "$aliases" --sample-metadata "$metadata" --outdir "$OUTDIR/common_input"
  python3 analysis_v2/scripts/shard_biomarker_model_input.py --profile-manifest "$OUTDIR/common_input/biomarker_profile_manifest.tsv" \
    --abundance-long "$OUTDIR/common_input/biomarker_abundance_long.tsv" --outdir "$inputs"
  ;;
 paired_map)
  i=$(printf '%04d' "$SLURM_ARRAY_TASK_ID"); d="$inputs/paired/$i"
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_paired_biomarker_models.R \
    --profile-manifest "$d/profile_manifest.tsv" --abundance-long "$d/abundance_long.tsv" --outdir "$d/pooled" --condition-mode pooled
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_paired_biomarker_models.R \
    --profile-manifest "$d/profile_manifest.tsv" --abundance-long "$d/abundance_long.tsv" --outdir "$d/stratified" --condition-mode stratified
  ;;
 disease_map)
  i=$(printf '%04d' "$SLURM_ARRAY_TASK_ID"); d="$inputs/disease/$i"
  apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_disease_biomarker_models.R \
    --profile-manifest "$d/profile_manifest.tsv" --abundance-long "$d/abundance_long.tsv" --outdir "$d/model"
  ;;
 paired_reduce)
  for mode in pooled stratified; do
    dest="$OUTDIR/models/artificial_${mode}"; mkdir -p "$dest"
    python3 analysis_v2/scripts/merge_biomarker_model_shards.py --shard-manifest "$inputs/shard_manifest.tsv" --kind paired --model-subdir "$mode" --outdir "$dest/models"
    python3 analysis_v2/scripts/evaluate_biomarker_propagation.py --calls "$dest/models/paired_da_results.tsv" --aliases "$aliases" --spike-panel "$panel" --outdir "$dest/evaluation" --q-thresholds 0.05,0.10
    printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$dest/DEVELOPMENT_ONLY.txt"
    printf 'analysis\tpaired_biomarker_propagation\nanalysis_status\tDEVELOPMENT_ONLY\nstatus\tPASS\n' > "$dest/SUCCESS"
  done
  ;;
 disease_reduce)
  dest="$OUTDIR/models/disease"; mkdir -p "$dest"
  python3 analysis_v2/scripts/merge_biomarker_model_shards.py --shard-manifest "$inputs/shard_manifest.tsv" --kind disease --model-subdir model --outdir "$dest/models"
  python3 analysis_v2/scripts/evaluate_disease_biomarker_propagation.py --calls "$dest/models/primary_disease_da_results.tsv" --aliases "$aliases" --spike-panel "$panel" --outdir "$dest/evaluation" --q-thresholds 0.05,0.10
  printf 'status\tDEVELOPMENT_ONLY\nuse_for_manuscript\tNO\n' > "$dest/DEVELOPMENT_ONLY.txt"
  printf 'analysis\tnative_disease_biomarker_propagation\nanalysis_status\tDEVELOPMENT_ONLY\nstatus\tPASS\n' > "$dest/SUCCESS"
  ;;
 artificial_report) env PAIRED_RUN="$OUTDIR/models/artificial_pooled" SECONDARY_PAIRED_RUN="$OUTDIR/models/artificial_stratified" ANALYSIS_SIF="$ANALYSIS_SIF" REPORT_STATUS=DEVELOPMENT_ONLY OUTDIR="$OUTDIR/reports/artificial" bash analysis_v2/run_artificial_biomarker_report.sh ;;
 disease_report) env DISEASE_RUN="$OUTDIR/models/disease" ANALYSIS_SIF="$ANALYSIS_SIF" REPORT_STATUS=DEVELOPMENT_ONLY OUTDIR="$OUTDIR/reports/disease" bash analysis_v2/run_disease_biomarker_report.sh ;;
 linkage) env ENDPOINTS_FILE="$DEV_INPUT_ROOT/endpoints/paired_endpoints.tsv" ENDPOINTS_SUCCESS="$DEV_INPUT_ROOT/endpoints/SUCCESS" PAIRED_RUN="$OUTDIR/models/artificial_pooled" SECONDARY_PAIRED_RUN="$OUTDIR/models/artificial_stratified" ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEVELOPMENT_ONLY OUTDIR="$OUTDIR/reports/calibration_linkage" bash analysis_v2/run_calibration_biomarker_linkage.sh ;;
 seal)
  req=(models/artificial_pooled models/artificial_stratified models/disease reports/artificial reports/disease reports/calibration_linkage)
  successes=(); for d in "${req[@]}"; do test -s "$OUTDIR/$d/SUCCESS" || { echo "[ERROR] missing $d/SUCCESS" >&2; exit 1; }; successes+=("$OUTDIR/$d/SUCCESS"); done
  sha256sum "$canonical" "$DEV_INPUT_ROOT/endpoints/paired_endpoints.tsv" "$metadata" "${successes[@]}" > "$OUTDIR/provenance/development_run.sha256"
  printf 'analysis\tlegacy_feng_zeller_biomarker_mapreduce\nstatus\tDEVELOPMENT_ONLY\ncompletion\tPASS\n' > "$OUTDIR/SUCCESS" ;;
 *) echo "[ERROR] unknown stage $PIPELINE_STAGE" >&2; exit 2;;
esac
