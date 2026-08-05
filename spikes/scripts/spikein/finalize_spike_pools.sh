#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<'EOF'
Usage:
  finalize_spike_pools.sh --env spikein.env [--delete-raw]

Validates both final pool mates for every panel target, writes an inventory and
SHA-256 manifest, verifies that manifest, and optionally deletes only the
redundant <label>.raw_1.fq and <label>.raw_2.fq files.

Raw deletion occurs only after all final pools pass validation and checksum
verification. Temporary/interrupted ART files are treated as an error.
EOF
}

ENV=""
DELETE_RAW=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV="$2"; shift 2 ;;
    --delete-raw) DELETE_RAW=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$ENV" && -f "$ENV" ]] || { usage; exit 2; }
# shellcheck disable=SC1090
source "$ENV"
: "${SPIKE_PANEL:?set SPIKE_PANEL in env}"
: "${POOLS_DIR:?set POOLS_DIR in env}"
[[ -s "$SPIKE_PANEL" ]] || { echo "[ERROR] Missing/empty panel: $SPIKE_PANEL" >&2; exit 1; }
[[ -d "$POOLS_DIR" ]] || { echo "[ERROR] Missing pool directory: $POOLS_DIR" >&2; exit 1; }

if compgen -G "$POOLS_DIR/*.art*.fq" >/dev/null; then
  echo "[ERROR] Interrupted ART temporary files remain in $POOLS_DIR" >&2
  printf '  %s\n' "$POOLS_DIR"/*.art*.fq >&2
  exit 1
fi

declare -a pool_files=()
declare -a raw_files=()
while IFS=$'\t' read -r label _taxon _assembly _fasta _weight _url; do
  [[ "$label" == "label" || -z "$label" ]] && continue
  for mate in 1 2; do
    pool="$POOLS_DIR/${label}.pool_${mate}.fq"
    raw="$POOLS_DIR/${label}.raw_${mate}.fq"
    [[ -s "$pool" ]] || { echo "[ERROR] Missing/empty final pool: $pool" >&2; exit 1; }
    pool_files+=("$(basename "$pool")")
    [[ -e "$raw" ]] && raw_files+=("$raw")
  done
done < "$SPIKE_PANEL"

(( ${#pool_files[@]} > 0 && ${#pool_files[@]} % 2 == 0 )) || {
  echo "[ERROR] No complete final pool set was discovered" >&2
  exit 1
}

inventory_tmp="$POOLS_DIR/pool_inventory.tsv.tmp.$$"
checksums_tmp="$POOLS_DIR/pool_files.sha256.tmp.$$"
trap 'rm -f "$inventory_tmp" "$checksums_tmp"' EXIT

{
  printf 'file\tbytes\n'
  for file in "${pool_files[@]}"; do
    printf '%s\t%s\n' "$file" "$(stat -c '%s' "$POOLS_DIR/$file")"
  done
} > "$inventory_tmp"

(
  cd "$POOLS_DIR"
  sha256sum "${pool_files[@]}"
) > "$checksums_tmp"

(
  cd "$POOLS_DIR"
  sha256sum -c "$(basename "$checksums_tmp")"
)

mv -f "$inventory_tmp" "$POOLS_DIR/pool_inventory.tsv"
mv -f "$checksums_tmp" "$POOLS_DIR/pool_files.sha256"
trap - EXIT

echo "[PASS] Validated and checksummed ${#pool_files[@]} final pool files"
echo "[OK] Inventory: $POOLS_DIR/pool_inventory.tsv"
echo "[OK] Checksums: $POOLS_DIR/pool_files.sha256"

if (( DELETE_RAW == 1 )); then
  if (( ${#raw_files[@]} > 0 )); then
    raw_bytes=0
    for raw in "${raw_files[@]}"; do
      raw_bytes=$((raw_bytes + $(stat -c '%s' "$raw")))
    done
    rm -f -- "${raw_files[@]}"
    echo "[OK] Deleted ${#raw_files[@]} verified-redundant raw ART files ($raw_bytes bytes)"
  else
    echo "[SKIP] No raw ART files remain"
  fi
else
  echo "[KEEP] Raw ART files retained; rerun with --delete-raw after reviewing the checksum result."
fi
