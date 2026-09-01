#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT"

: "${UHGG_METADATA:?Set UHGG_METADATA to genomes-all_metadata.tsv}"
: "${UHGG_MASH_SKETCH:?Set UHGG_MASH_SKETCH to all_genomes.msh}"
: "${METAPHLAN_PKL:?Set METAPHLAN_PKL to the frozen MetaPhlAn database pickle}"
: "${NCBI_REPORT_DIR:?Set NCBI_REPORT_DIR to the production assembly provenance directory}"
: "${PHASE3_OUT:?Set PHASE3_OUT to an empty/new output directory}"

SPIKE_PANEL="${SPIKE_PANEL:-$ROOT/spikes/spike_panel.tsv}"
ALIASES="${ALIASES:-$ROOT/examples/spike_taxon_aliases.csv}"
MASH_BIN="${MASH_BIN:-$ROOT/mash_in_container.sh}"

for path in "$SPIKE_PANEL" "$ALIASES" "$UHGG_METADATA" "$UHGG_MASH_SKETCH" "$METAPHLAN_PKL"; do
  [[ -s "$path" ]] || { echo "[ERROR] Missing or empty input: $path" >&2; exit 1; }
done
[[ ! -e "$PHASE3_OUT/SUCCESS" ]] || {
  echo "[ERROR] Refusing to overwrite a completed audit: $PHASE3_OUT" >&2
  exit 1
}

mkdir -p "$PHASE3_OUT/assembly_integrity" "$PHASE3_OUT/database_representation"
python3 scripts/audit_target_assemblies.py \
  --spike-panel "$SPIKE_PANEL" \
  --ncbi-report-dir "$NCBI_REPORT_DIR" \
  --outdir "$PHASE3_OUT/assembly_integrity"
python3 scripts/build_reference_representation_table.py \
  --spike-panel "$SPIKE_PANEL" \
  --aliases "$ALIASES" \
  --uhgg-metadata "$UHGG_METADATA" \
  --uhgg-mash-sketch "$UHGG_MASH_SKETCH" \
  --metaphlan-pkl "$METAPHLAN_PKL" \
  --mash-bin "$MASH_BIN" \
  --outdir "$PHASE3_OUT/database_representation"

find "$PHASE3_OUT" -type f ! -name 'PHASE3_OUTPUTS.sha256' ! -name SUCCESS \
  -print0 | sort -z | xargs -0 sha256sum > "$PHASE3_OUT/PHASE3_OUTPUTS.sha256"
touch "$PHASE3_OUT/SUCCESS"
echo "[PASS] Phase 3 target validation completed: $PHASE3_OUT"
