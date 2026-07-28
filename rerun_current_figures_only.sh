#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-$PWD}"
: "${RUN_ROOT:?Set RUN_ROOT to the completed analysis directory}"
R_BIN="${R_BIN:-$PROJECT/rscript_in_container.sh}"

cd "$PROJECT"
PROJECT="$(pwd -P)"

if [[ "$RUN_ROOT" != /* ]]; then
  RUN_ROOT="$PROJECT/$RUN_ROOT"
fi
[[ -d "$RUN_ROOT" ]] || {
  echo "[ERROR] Completed RUN_ROOT does not exist: $RUN_ROOT" >&2
  exit 1
}
[[ -x "$R_BIN" ]] || {
  echo "[ERROR] R_BIN is not executable: $R_BIN" >&2
  exit 1
}

required=(
  "$RUN_ROOT/manuscript_figures/current/source_panels"
  "$PROJECT/scripts/99_make_current_manuscript_figures.R"
)
for path in "${required[@]}"; do
  [[ -e "$path" ]] || {
    echo "[ERROR] Missing required path: $path" >&2
    exit 1
  }
done

backup="$RUN_ROOT/manuscript_figures/current/blank_composites_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$backup"
for stem in \
  Fig2_baseline_profiler_discordance \
  Supplementary_Fig_B2_baseline_all_10_targets \
  Fig3_independent_spike_recovery \
  Fig4_biomarker_recoverability \
  Fig5_community_spike_recovery_and_DA
do
  for ext in pdf png; do
    path="$RUN_ROOT/manuscript_figures/current/$stem.$ext"
    [[ -f "$path" ]] && mv "$path" "$backup/"
  done
done

echo "[INFO] Project:             $PROJECT"
echo "[INFO] Completed run:       $RUN_ROOT"
echo "[INFO] Existing composites: $backup"
echo "[INFO] Reusing source panels; analyses will not be rerun."

"$R_BIN" "$PROJECT/scripts/99_make_current_manuscript_figures.R" \
  --project-root "$PROJECT" \
  --run-dir "$RUN_ROOT" \
  --metadata metadata_w_study.tsv \
  --kraken kraken2_bracken_merged_unspecified.csv \
  --metaphlan metaphlan4_merged_unspecified.csv \
  --rerun-panels false

echo "[PASS] Current manuscript composites regenerated from existing source panels."
echo "[INFO] Output: $RUN_ROOT/manuscript_figures/current"
