#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"; cd "$ROOT"
: "${LEGACY_INPUT_ROOT:?Set the validated legacy Feng/Zeller input root}"
: "${YACHIDA_ENV:?Set the strict-production Yachida environment}"
: "${YACHIDA_MANIFEST:?Set the complete Yachida manifest}"
: "${YACHIDA_INDEPENDENT_MANIFEST:?Set the Yachida independent subset manifest}"
: "${OUTDIR:?Set a new combined development input directory}"
source "$YACHIDA_ENV"
test -s "$YACHIDA_STATE_DIR/production_seal/SUCCESS" || { echo '[ERROR] Yachida production is not sealed' >&2; exit 1; }
test ! -e "$OUTDIR" || { echo "[ERROR] OUTDIR exists: $OUTDIR" >&2; exit 1; }
mkdir -p "$OUTDIR/yachida"
python3 analysis_v2/scripts/build_crc_cohort_canonical_input.py \
  --cohort yachida --manifest "$YACHIDA_MANIFEST" \
  --independent-manifest "$YACHIDA_INDEPENDENT_MANIFEST" \
  --results-root "$PERSISTENT_RESULTS_ROOT" --spike-panel spikes/spike_panel.tsv \
  --aliases examples/spike_taxon_aliases.csv --outdir "$OUTDIR/yachida"
python3 analysis_v2/scripts/combine_canonical_development_inputs.py \
  --input "$LEGACY_INPUT_ROOT/canonical_input.tsv" \
  --input "$OUTDIR/yachida/canonical_input.tsv" --outdir "$OUTDIR/combined"
python3 analysis_v2/scripts/derive_paired_endpoints.py \
  --input "$OUTDIR/combined/canonical_input.tsv" --outdir "$OUTDIR/combined/endpoints"
cp "$OUTDIR/combined/DEVELOPMENT_ONLY.txt" "$OUTDIR/DEVELOPMENT_ONLY.txt"
printf 'combined_input\t%s\nyachida_seal\t%s\n' \
  "$OUTDIR/combined/canonical_input.tsv" "$YACHIDA_STATE_DIR/production_seal/SUCCESS" > "$OUTDIR/SUCCESS"
echo "[PASS] Three-cohort development input ready: $OUTDIR/combined"
