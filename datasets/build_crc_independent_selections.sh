#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
OUTPUT_ROOT="${1:-$ROOT/work/selection_preview}"
SELECTOR="$ROOT/scripts/select_samples_deterministically.py"
AUDITOR="$ROOT/datasets/audit_independent_selection.py"
ELIGIBILITY_BUILDER="$ROOT/datasets/build_crc_eligible_manifest.py"

build_one() {
  local cohort="$1"
  local source="$2"
  local study_name="$3"
  local seed="ground-truth-taxonomic-calibration-${cohort}-independent-v1"
  local destination="$OUTPUT_ROOT/$cohort"
  local eligible="$destination/eligible_cohort.tsv"
  local selection="$destination/independent_10_per_condition.tsv"

  mkdir -p "$destination"
  python3 "$ELIGIBILITY_BUILDER" \
    --metadata "$source" \
    --output "$eligible" \
    --study-name "$study_name"
  python3 "$SELECTOR" \
    --manifest "$eligible" \
    --output "$selection" \
    --per-condition 10 \
    --selection-seed "$seed" \
    --id-column Name \
    --condition-column "Study condition"
  python3 "$AUDITOR" \
    --eligible-manifest "$eligible" \
    --selection "$selection" \
    --output-dir "$destination/audit"
}

build_one fengq "$ROOT/datasets/fengq/Public_study__FengQ_2015.tsv" "Public_study__FengQ_2015"
build_one zellerg "$ROOT/datasets/zellerg/Public_study__ZellerG_2014.tsv" "Public_study__ZellerG_2014"

echo "[PASS] Built deterministic cohort selections under: $OUTPUT_ROOT"
