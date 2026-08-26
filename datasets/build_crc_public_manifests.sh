#!/usr/bin/env bash
# Rebuild the frozen Feng and Zeller public production manifests.
set -euo pipefail
IFS=$'\n\t'

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
RETRIEVED_UTC=2026-08-26
cd "$ROOT"

build_one() {
  local cohort="$1" metadata="$2" study_name="$3" inventory="$4"
  local ena_accession="$5" seed="$6"
  local output="datasets/$cohort/manifests"
  local eligible="$output/eligible_cohort.tsv"
  local selection="$output/independent_10_per_condition.tsv"
  local ena_report="$output/${ena_accession}.read_run.tsv"

  [[ -s "$ena_report" ]] || {
    echo "[ERROR] Missing frozen ENA snapshot: $ena_report" >&2
    exit 1
  }
  mkdir -p "$output/selection_audit"

  python3 datasets/build_crc_eligible_manifest.py \
    --metadata "$metadata" \
    --output "$eligible" \
    --study-name "$study_name"

  python3 scripts/select_samples_deterministically.py \
    --manifest "$eligible" \
    --output "$selection" \
    --per-condition 10 \
    --selection-seed "$seed" \
    --id-column Name \
    --condition-column "Study condition"

  python3 datasets/audit_independent_selection.py \
    --eligible-manifest "$eligible" \
    --selection "$selection" \
    --output-dir "$output/selection_audit"

  python3 datasets/build_crc_production_manifest.py \
    --eligible-manifest "$eligible" \
    --url-inventory "$inventory" \
    --ena-report "$ena_report" \
    --independent-selection "$selection" \
    --ena-study-accession "$ena_accession" \
    --ena-report-retrieved-date "$RETRIEVED_UTC" \
    --output "$output/production_manifest.tsv"

  (
    cd "$output"
    sha256sum "${ena_accession}.read_run.tsv" > "${ena_accession}.read_run.tsv.sha256"
  )
}

build_one \
  fengq \
  datasets/fengq/Public_study__FengQ_2015.tsv \
  Public_study__FengQ_2015 \
  datasets/fengq/fengq_wgets.sh \
  PRJEB7774 \
  ground-truth-taxonomic-calibration-fengq-independent-v1

build_one \
  zellerg \
  datasets/zellerg/Public_study__ZellerG_2014.tsv \
  Public_study__ZellerG_2014 \
  datasets/zellerg/zellerg_wgets.sh \
  PRJEB6070 \
  ground-truth-taxonomic-calibration-zellerg-independent-v1

echo "[PASS] Rebuilt frozen public Feng and Zeller manifests"
