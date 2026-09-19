#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)";cd "$ROOT"
: "${PROFILE_MANIFEST:?Set PROFILE_MANIFEST}"
: "${ABUNDANCE_LONG:?Set ABUNDANCE_LONG}"
: "${PAIRED_ENDPOINTS:?Set PAIRED_ENDPOINTS}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF}"
: "${OUTDIR:?Set OUTDIR to a new directory}"
: "${ANALYSIS_STATUS:?Set ANALYSIS_STATUS}"
if [[ "$ANALYSIS_STATUS" != DEVELOPMENT_ONLY ]];then echo "[ERROR] External rule is not frozen; only DEVELOPMENT_ONLY is allowed" >&2;exit 1;fi
for p in "$PROFILE_MANIFEST" "$ABUNDANCE_LONG" "$PAIRED_ENDPOINTS" "$ANALYSIS_SIF";do [[ -s "$p" ]]||{ echo "[ERROR] missing input: $p" >&2;exit 1;};done
if [[ -e "$OUTDIR" ]]&&[[ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]];then echo "[ERROR] OUTDIR must be new or empty" >&2;exit 1;fi
python3 analysis_v2/tests/test_cross_cohort_abundance_calibration.py
python3 analysis_v2/tests/test_calibration_restoration.py
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_disease_biomarker_models.R
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/tests/test_calibration_restoration_figures.R
python3 analysis_v2/scripts/validate_analysis_policy.py --policy analysis_v2/ANALYSIS_POLICY.tsv --outdir "$OUTDIR/analysis_policy"
python3 analysis_v2/scripts/calibrate_abundance_cross_cohort.py \
 --profile-manifest "$PROFILE_MANIFEST" --abundance-long "$ABUNDANCE_LONG" \
 --paired-endpoints "$PAIRED_ENDPOINTS" --outdir "$OUTDIR/calibration"
while IFS=$'\t' read -r transfer training validation directory active profiles;do
 [[ "$transfer" == transfer ]]&&continue
 mkdir -p "$directory/models_uncorrected" "$directory/models_corrected"
 apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_disease_biomarker_models.R \
  --profile-manifest "$directory/uncorrected_profile_manifest.tsv" --abundance-long "$directory/uncorrected_abundance_long.tsv" --outdir "$directory/models_uncorrected"
 apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/fit_disease_biomarker_models.R \
  --profile-manifest "$directory/corrected_profile_manifest.tsv" --abundance-long "$directory/corrected_abundance_long.tsv" --outdir "$directory/models_corrected"
 python3 analysis_v2/scripts/evaluate_calibration_restoration.py \
  --uncorrected-results "$directory/models_uncorrected/primary_disease_da_results.tsv" \
  --corrected-results "$directory/models_corrected/primary_disease_da_results.tsv" \
  --corrected-manifest "$directory/corrected_profile_manifest.tsv" --outdir "$directory/restoration"
done < "$OUTDIR/calibration/transfer_ledger.tsv"
apptainer exec --cleanenv --pwd "$ROOT" "$ANALYSIS_SIF" Rscript analysis_v2/scripts/plot_calibration_restoration.R \
 --calibration-root "$OUTDIR/calibration" --outdir "$OUTDIR/report"
for p in "$OUTDIR/calibration/SUCCESS" "$OUTDIR/report/SUCCESS" "$OUTDIR/report/figures/distortion_coefficient_atlas.pdf" "$OUTDIR/report/figures/abundance_error_restoration.pdf" "$OUTDIR/report/figures/biomarker_call_restoration.pdf" "$OUTDIR/report/figures/disease_effect_restoration.pdf";do [[ -s "$p" ]]||{ echo "[ERROR] missing output: $p" >&2;exit 1;};done
printf 'status=DEVELOPMENT_ONLY\nautomatic_feature_removal=NO\n' > "$OUTDIR/DEVELOPMENT_ONLY.txt";printf 'status=PASS\n' > "$OUTDIR/SUCCESS"
echo "[PASS] Cross-cohort abundance calibration and biomarker restoration: $OUTDIR"
