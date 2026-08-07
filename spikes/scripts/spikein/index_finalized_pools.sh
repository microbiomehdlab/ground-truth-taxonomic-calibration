#!/usr/bin/env bash
# Count finalized paired pools once so streaming tasks never rescan them for size.
set -euo pipefail
IFS=$'\n\t'

[[ ${1:-} == --env && -n ${2:-} && -f $2 ]] || {
  echo "Usage: $0 --env spikein.env" >&2
  exit 2
}
source "$2"
: "${SPIKE_PANEL:?set SPIKE_PANEL}"
: "${POOLS_DIR:?set POOLS_DIR}"

[[ -s "$POOLS_DIR/pool_files.sha256" ]] || { echo "[ERROR] Final pool checksums are missing" >&2; exit 1; }
echo "[INFO] Using the checksum manifest emitted by finalize_spike_pools.sh"
echo "[INFO] Pool structure and synchronized pair counts will be checked in one scan"

tmp="$POOLS_DIR/pool_pair_counts.tsv.tmp.$$"
trap 'rm -f "$tmp"' EXIT
printf 'label\tmate1_pairs\tmate2_pairs\tpool_pairs\n' > "$tmp"
while IFS=$'\t' read -r label _taxon _assembly _fasta _weight _url; do
  [[ "$label" == label || -z "$label" ]] && continue
  p1="$POOLS_DIR/${label}.pool_1.fq"
  p2="$POOLS_DIR/${label}.pool_2.fq"
  n1="$(awk 'END {print int(NR/4)}' "$p1")"
  n2="$(awk 'END {print int(NR/4)}' "$p2")"
  [[ "$n1" -gt 0 && "$n1" == "$n2" ]] || {
    echo "[ERROR] Invalid paired pool size for $label: R1=$n1 R2=$n2" >&2
    exit 1
  }
  printf '%s\t%s\t%s\t%s\n' "$label" "$n1" "$n2" "$n1" >> "$tmp"
  echo "[OK] $label: $n1 pairs"
done < "$SPIKE_PANEL"
mv -f "$tmp" "$POOLS_DIR/pool_pair_counts.tsv"
sha256sum "$POOLS_DIR/pool_pair_counts.tsv" > "$POOLS_DIR/pool_pair_counts.tsv.sha256"
sha256sum "$POOLS_DIR/pool_files.sha256" > "$POOLS_DIR/pool_files.sha256.sha256"
trap - EXIT
echo "[PASS] Finalized pool pair-count index written: $POOLS_DIR/pool_pair_counts.tsv"
