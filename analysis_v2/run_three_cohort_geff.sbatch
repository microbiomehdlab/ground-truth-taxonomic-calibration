#!/bin/bash
#SBATCH --job-name=three_cohort_geff
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=three_cohort_geff.%j.out
#SBATCH --error=three_cohort_geff.%j.err
set -euo pipefail
: "${PROJECT:?Export repository root}"
: "${GEFF_ROOT:?Export a fresh output root}"
: "${GENOME_SIZE_MAPPING:?Export the audited production-database length mapping}"
: "${EXPECTED_MAPPING_SHA256:?Export its previously verified SHA256}"
: "${ANALYSIS_SIF:?Export analysis image}"
cd "$PROJECT"
test ! -e "$GEFF_ROOT"
test -s "$ANALYSIS_SIF"
test -s "$GENOME_SIZE_MAPPING"
test "$(sha256sum "$GENOME_SIZE_MAPPING" | cut -d' ' -f1)" = "$EXPECTED_MAPPING_SHA256"
sha256sum -c --quiet analysis_v2/taxon_identity_freeze.sha256
appt=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF")
for cohort in yachida feng zeller; do
  (
    source "config/$cohort.strict-production.env"
    case "$cohort" in
      yachida) state="$YACHIDA_STATE_DIR"; expected=201 ;;
      feng) state="$CRC_STATE_DIR"; expected=154 ;;
      zeller) state="$CRC_STATE_DIR"; expected=156 ;;
    esac
    "${appt[@]}" python3 -B analysis_v2/scripts/prepare_sealed_metaphlan_baselines.py \
      --cohort "$cohort" --state-dir "$state" --results-root "$PERSISTENT_RESULTS_ROOT" \
      --expected-samples "$expected" --outdir "$GEFF_ROOT/$cohort/selection"
    "${appt[@]}" python3 -B analysis_v2/scripts/compute_effective_genome_size.py \
      --manifest "$GEFF_ROOT/$cohort/selection/baseline_manifest.tsv" \
      --genome-sizes "$GENOME_SIZE_MAPPING" --profile-rank sgb --match-on lineage \
      --abundance-unit percent --min-coverage 0.95 --max-excluded-samples 0 \
      --outdir "$GEFF_ROOT/$cohort/geff_primary_cov95"
    test -s "$GEFF_ROOT/$cohort/geff_primary_cov95/SUCCESS"
    (cd "$GEFF_ROOT/$cohort/geff_primary_cov95"; sha256sum -c --quiet effective_genome_size.sha256)
  )
done
"${appt[@]}" python3 -B analysis_v2/scripts/build_target_genome_sizes.py \
  --spike-panel spikes/spike_panel.tsv --fasta-root "$PROJECT" --outdir "$GEFF_ROOT/targets"
git rev-parse HEAD > "$GEFF_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$GEFF_ROOT/analysis_image.sha256"
sha256sum "$GENOME_SIZE_MAPPING" > "$GEFF_ROOT/mapping.sha256"
sha256sum analysis_v2/scripts/prepare_sealed_metaphlan_baselines.py \
  analysis_v2/scripts/compute_effective_genome_size.py \
  analysis_v2/scripts/build_target_genome_sizes.py analysis_v2/run_three_cohort_geff.sbatch \
  > "$GEFF_ROOT/code.sha256"
printf 'status\tPASS\nsamples\t511\n' > "$GEFF_ROOT/SUCCESS"
(cd "$GEFF_ROOT"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Sealed-baseline G_eff for all three cohorts: $GEFF_ROOT"
