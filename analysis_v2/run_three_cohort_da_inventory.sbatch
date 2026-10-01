#!/bin/bash
#SBATCH --job-name=da_inventory
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=da_inventory.%j.out
#SBATCH --error=da_inventory.%j.err
set -euo pipefail
: "${PROJECT:?Export repository root}"
: "${ANALYSIS_SIF:?Export pinned analysis image}"
: "${DA_INVENTORY_ROOT:?Export fresh inventory root}"
: "${PREFLIGHT_STAMP:?Export definitive shared preflight timestamp}"
cd "$PROJECT"
test ! -e "$DA_INVENTORY_ROOT"
test -s "$ANALYSIS_SIF"
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
    canonical="$PROJECT/work/analysis_v2_${cohort}_shared_preflight_${PREFLIGHT_STAMP}/canonical/canonical_input.tsv"
    "${appt[@]}" python3 -B analysis_v2/scripts/build_da_design_inventory.py \
      --canonical "$canonical" --state-dir "$state" --cohort "$cohort" \
      --expected-samples "$expected" --repetitions 1000 \
      --outdir "$DA_INVENTORY_ROOT/$cohort"
    test -s "$DA_INVENTORY_ROOT/$cohort/SUCCESS"
    (cd "$DA_INVENTORY_ROOT/$cohort"; sha256sum -c --quiet SHA256SUMS)
    cat "$DA_INVENTORY_ROOT/$cohort/eligibility.tsv"
  )
done
git rev-parse HEAD > "$DA_INVENTORY_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$DA_INVENTORY_ROOT/analysis_image.sha256"
sha256sum analysis_v2/scripts/build_da_design_inventory.py analysis_v2/scripts/da_allocations.py \
  analysis_v2/scripts/build_biomarker_abundance_input.py analysis_v2/scripts/audit_bracken_denominators.py \
  analysis_v2/run_three_cohort_da_inventory.sbatch > "$DA_INVENTORY_ROOT/code.sha256"
"${appt[@]}" python3 --version > "$DA_INVENTORY_ROOT/python_version.txt"
printf 'stage\tDA2_DA3_INPUT_AND_ALLOCATION_INVENTORY\nfull_pool_draws\t1000\nsubset_assignments\t252\nDA1_eligibility\tPENDING\nmodels_run\t0\n' > "$DA_INVENTORY_ROOT/design_stage.tsv"
printf 'status\tPASS_DESIGN_INVENTORY\ndefinitive_model_readiness\tNOT_ESTABLISHED\n' > "$DA_INVENTORY_ROOT/SUCCESS"
(cd "$DA_INVENTORY_ROOT"; find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS; sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Three-cohort DA inventory; no models run: $DA_INVENTORY_ROOT"
