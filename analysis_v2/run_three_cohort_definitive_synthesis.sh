#!/usr/bin/env bash
# Fail-closed bridge from three definitive cohort packages to meta-analysis.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${YACHIDA_RUN:?Set YACHIDA_RUN to its definitive package}"
: "${FENG_RUN:?Set FENG_RUN to its definitive package}"
: "${ZELLER_RUN:?Set ZELLER_RUN to its definitive package}"
: "${ANALYSIS_SIF:?Set ANALYSIS_SIF}"
: "${OUTDIR:?Set OUTDIR to a new synthesis directory}"
for item in "yachida:$YACHIDA_RUN" "feng:$FENG_RUN" "zeller:$ZELLER_RUN"; do
  cohort="${item%%:*}"; run="${item#*:}"
  for path in "$run/SUCCESS" "$run/readiness/SUCCESS" "$run/models/disease/SUCCESS" \
    "$run/models/disease/models/primary_disease_da_results.tsv" "$run/provenance/definitive_run.sha256"; do
    [[ -s "$path" ]] || { echo "[ERROR] Incomplete definitive $cohort package: $path" >&2; exit 1; }
  done
  [[ ! -e "$run/DEVELOPMENT_ONLY.txt" ]] || { echo "[ERROR] $cohort package is development-only" >&2; exit 1; }
done
COHORT_RESULT_FILES="$YACHIDA_RUN/models/disease/models/primary_disease_da_results.tsv:$FENG_RUN/models/disease/models/primary_disease_da_results.tsv:$ZELLER_RUN/models/disease/models/primary_disease_da_results.tsv" \
ANALYSIS_SIF="$ANALYSIS_SIF" ANALYSIS_STATUS=DEFINITIVE OUTDIR="$OUTDIR" \
  bash analysis_v2/run_cross_cohort_synthesis.sh
sha256sum "$YACHIDA_RUN/SUCCESS" "$FENG_RUN/SUCCESS" "$ZELLER_RUN/SUCCESS" \
  "$OUTDIR/SUCCESS" >> "$OUTDIR/provenance/run_inputs_and_primary_outputs.sha256"
echo "[PASS] Three-cohort definitive synthesis sealed: $OUTDIR"
