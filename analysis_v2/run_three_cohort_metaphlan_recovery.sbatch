#!/bin/bash
#SBATCH --job-name=metaphlan_recovery
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=metaphlan_recovery.%j.out
#SBATCH --error=metaphlan_recovery.%j.err
set -euo pipefail
: "${PROJECT:?Export project root}"
: "${GEFF_ROOT:?Export verified 511-sample reference bundle}"
: "${METAPHLAN_RECOVERY_ROOT:?Export a fresh output root}"
: "${ANALYSIS_SIF:?Export analysis image}"
: "${PREFLIGHT_STAMP:?Export definitive canonical preflight timestamp}"
cd "$PROJECT"
test ! -e "$METAPHLAN_RECOVERY_ROOT"
test -s "$GEFF_ROOT/SUCCESS"
(cd "$GEFF_ROOT"; sha256sum -c --quiet SHA256SUMS)
sha256sum -c --quiet analysis_v2/taxon_identity_freeze.sha256
appt=(apptainer exec --cleanenv --bind /mnt/nfs:/mnt/nfs --pwd "$PROJECT" "$ANALYSIS_SIF")
for cohort in yachida feng zeller; do
  case "$cohort" in yachida) expected=201 ;; feng) expected=154 ;; zeller) expected=156 ;; esac
  canonical="$PROJECT/work/analysis_v2_${cohort}_shared_preflight_${PREFLIGHT_STAMP}/canonical/canonical_input.tsv"
  selection="$GEFF_ROOT/$cohort/selection"
  geff="$GEFF_ROOT/$cohort/geff_primary_cov95/effective_genome_size.tsv"
  targets="$GEFF_ROOT/targets/target_genome_sizes.tsv"
  "${appt[@]}" python3 -B analysis_v2/scripts/check_geff_canonical_baselines.py \
    --canonical "$canonical" --selection "$selection" --cohort "$cohort" --expected-samples "$expected"
  "${appt[@]}" python3 -B analysis_v2/scripts/validate_metaphlan_reference_inputs.py \
    --canonical "$canonical" --targets "$targets" --effective "$geff"
  # A single primary derivation preserves unchanged read-reference fields too.
  # It emits BOTH profiler rows; Bracken here remains the legacy comparator,
  # not the separately established all-input Bracken primary measurement.
  "${appt[@]}" python3 -B analysis_v2/scripts/derive_paired_endpoints.py \
    --input "$canonical" --outdir "$METAPHLAN_RECOVERY_ROOT/$cohort/endpoints" \
    --metaphlan-reference genome_equivalent --target-genome-sizes "$targets" \
    --effective-genome-sizes "$geff"
  test -s "$METAPHLAN_RECOVERY_ROOT/$cohort/endpoints/SUCCESS"
  (cd "$METAPHLAN_RECOVERY_ROOT/$cohort/endpoints"; sha256sum -c --quiet endpoint_derivation.sha256)
done
git rev-parse HEAD > "$METAPHLAN_RECOVERY_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$METAPHLAN_RECOVERY_ROOT/analysis_image.sha256"
sha256sum analysis_v2/scripts/check_geff_canonical_baselines.py \
  analysis_v2/scripts/validate_metaphlan_reference_inputs.py \
  analysis_v2/scripts/derive_paired_endpoints.py analysis_v2/run_three_cohort_metaphlan_recovery.sbatch \
  > "$METAPHLAN_RECOVERY_ROOT/code.sha256"
printf 'status\tPASS_ENDPOINTS\nreference\tgenome_equivalent_with_read_reference_fields\n' > "$METAPHLAN_RECOVERY_ROOT/SUCCESS"
(cd "$METAPHLAN_RECOVERY_ROOT"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Three-cohort recovery endpoints; no DA models run: $METAPHLAN_RECOVERY_ROOT"
