#!/usr/bin/env bash
# Build the upstream evidence package inside the pinned analysis container.
#
# It consumes the three verified production_seal_v2 directories and the sealed
# Yachida assembly-sensitivity experiment. It never mutates a source seal,
# never reruns profiling, and never copies raw reads, databases, container
# images, scratch files or scheduler logs.
#
# Every authoritative path is passed explicitly; nothing is discovered by
# searching work/.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT"
set -euo pipefail
IFS=$'\n\t'

: "${ANALYSIS_SIF:?export ANALYSIS_SIF to the pinned analysis image}"
: "${OUTDIR:?export OUTDIR to a new package directory}"
: "${YACHIDA_SEAL:?export YACHIDA_SEAL to its production_seal_v2 directory}"
: "${FENG_SEAL:?export FENG_SEAL to its production_seal_v2 directory}"
: "${ZELLER_SEAL:?export ZELLER_SEAL to its production_seal_v2 directory}"
: "${ASSEMBLY_SENSITIVITY_SEAL:?export ASSEMBLY_SENSITIVITY_SEAL to the experiment_seal directory}"
: "${YACHIDA_MANIFEST:?export YACHIDA_MANIFEST}"
: "${YACHIDA_INDEPENDENT_MANIFEST:?export YACHIDA_INDEPENDENT_MANIFEST}"
: "${FENG_MANIFEST:?export FENG_MANIFEST}"
: "${FENG_INDEPENDENT_MANIFEST:?export FENG_INDEPENDENT_MANIFEST}"
: "${ZELLER_MANIFEST:?export ZELLER_MANIFEST}"
: "${ZELLER_INDEPENDENT_MANIFEST:?export ZELLER_INDEPENDENT_MANIFEST}"
: "${PROVENANCE_METADATA:?export PROVENANCE_METADATA to upstream_provenance_metadata.tsv}"
: "${AUDIT_LEDGER:?export AUDIT_LEDGER to the manually supplied audit ledger}"

SPIKE_PANEL="${SPIKE_PANEL:-$ROOT/spikes/spike_panel.tsv}"
TAXON_ALIASES="${TAXON_ALIASES:-$ROOT/examples/spike_taxon_aliases.csv}"

test -s "$ANALYSIS_SIF" || { echo "[ERROR] missing image: $ANALYSIS_SIF" >&2; exit 1; }
[[ ! -e "$OUTDIR" || -z "$(find "$OUTDIR" -mindepth 1 -print -quit)" ]] || {
  echo "[ERROR] OUTDIR must be new or empty: $OUTDIR" >&2
  exit 1
}
for required in "$YACHIDA_SEAL" "$FENG_SEAL" "$ZELLER_SEAL" \
                "$ASSEMBLY_SENSITIVITY_SEAL"; do
  test -d "$required" || { echo "[ERROR] not a directory: $required" >&2; exit 1; }
  test -s "$required/SUCCESS" || {
    echo "[ERROR] source seal has no SUCCESS: $required" >&2; exit 1; }
  test ! -e "$required/AUDIT_IN_PROGRESS" || {
    echo "[ERROR] source seal is mid-audit: $required" >&2; exit 1; }
done
for required in "$YACHIDA_MANIFEST" "$YACHIDA_INDEPENDENT_MANIFEST" \
                "$FENG_MANIFEST" "$FENG_INDEPENDENT_MANIFEST" \
                "$ZELLER_MANIFEST" "$ZELLER_INDEPENDENT_MANIFEST" \
                "$SPIKE_PANEL" "$TAXON_ALIASES" "$PROVENANCE_METADATA" \
                "$AUDIT_LEDGER"; do
  test -s "$required" || { echo "[ERROR] missing or empty: $required" >&2; exit 1; }
done

COMMIT="$(git -C "$ROOT" rev-parse HEAD)"
IMAGE_SHA256="$(sha256sum "$ANALYSIS_SIF" | cut -d' ' -f1)"
echo "[INFO] analysis image sha256: $IMAGE_SHA256"
echo "[INFO] source commit: $COMMIT"

# Snapshot every source seal as a complete inventory, so a file that is
# added, removed or changed during the build is all detectable. Checking the
# original list alone would miss a newly added file.
SEAL_BEFORE="$(mktemp)"
SEAL_AFTER="$(mktemp)"
trap 'rm -f "$SEAL_BEFORE" "$SEAL_AFTER"' EXIT
seal_inventory() {
  find "$YACHIDA_SEAL" "$FENG_SEAL" "$ZELLER_SEAL" \
    "$ASSEMBLY_SENSITIVITY_SEAL" -type f -print0 \
    | LC_ALL=C sort -z | xargs -0 sha256sum
}
seal_inventory > "$SEAL_BEFORE"

apptainer exec --cleanenv --bind "$PWD:$PWD" --pwd "$PWD" "$ANALYSIS_SIF" \
  python3 analysis_v2/scripts/build_upstream_evidence_package.py \
    --yachida-seal "$YACHIDA_SEAL" \
    --yachida-manifest "$YACHIDA_MANIFEST" \
    --yachida-independent-manifest "$YACHIDA_INDEPENDENT_MANIFEST" \
    --feng-seal "$FENG_SEAL" \
    --feng-manifest "$FENG_MANIFEST" \
    --feng-independent-manifest "$FENG_INDEPENDENT_MANIFEST" \
    --zeller-seal "$ZELLER_SEAL" \
    --zeller-manifest "$ZELLER_MANIFEST" \
    --zeller-independent-manifest "$ZELLER_INDEPENDENT_MANIFEST" \
    --assembly-sensitivity-seal "$ASSEMBLY_SENSITIVITY_SEAL" \
    --spike-panel "$SPIKE_PANEL" \
    --taxon-aliases "$TAXON_ALIASES" \
    --provenance-metadata "$PROVENANCE_METADATA" \
    --audit-ledger "$AUDIT_LEDGER" \
    --source-commit "$COMMIT" \
    --analysis-image-sha256 "$IMAGE_SHA256" \
    --require-figure \
    --outdir "$OUTDIR"

# Compare complete before/after inventories: this detects a changed file, a
# removed file and a newly added one.
seal_inventory > "$SEAL_AFTER"
if ! diff -u "$SEAL_BEFORE" "$SEAL_AFTER" > /dev/null; then
  echo "[ERROR] a source seal changed during the build:" >&2
  diff -u "$SEAL_BEFORE" "$SEAL_AFTER" >&2 || true
  if [[ -d "$OUTDIR" ]]; then
    REJECTED_OUTDIR="${OUTDIR}.rejected-source-mutation"
    [[ ! -e "$REJECTED_OUTDIR" ]] || {
      echo "[ERROR] cannot quarantine rejected package; target exists: $REJECTED_OUTDIR" >&2
      exit 1
    }
    mv -- "$OUTDIR" "$REJECTED_OUTDIR"
    echo "[ERROR] rejected package quarantined at: $REJECTED_OUTDIR" >&2
  fi
  exit 1
fi
for seal in "$YACHIDA_SEAL" "$FENG_SEAL" "$ZELLER_SEAL" \
            "$ASSEMBLY_SENSITIVITY_SEAL"; do
  test ! -e "$seal/AUDIT_IN_PROGRESS" || {
    echo "[ERROR] a source seal became mid-audit during the build: $seal" >&2
    if [[ -d "$OUTDIR" ]]; then
      REJECTED_OUTDIR="${OUTDIR}.rejected-source-mutation"
      [[ ! -e "$REJECTED_OUTDIR" ]] && mv -- "$OUTDIR" "$REJECTED_OUTDIR"
    fi
    exit 1
  }
done

test -s "$OUTDIR/SUCCESS" || { echo "[ERROR] package has no SUCCESS" >&2; exit 1; }
test -s "$OUTDIR/SHA256SUMS" || { echo "[ERROR] package has no SHA256SUMS" >&2; exit 1; }
(cd "$OUTDIR" && sha256sum -c --quiet SHA256SUMS)
echo "[PASS] Upstream evidence package verified: $OUTDIR"
