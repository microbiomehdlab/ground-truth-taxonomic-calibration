#!/usr/bin/env bash
set -euo pipefail

# Download only the UHGG v2.0.2 companion assets required for Supplementary
# Table A7. The Kraken2 index and MetaPhlAn database are not downloaded here.

CACHE_DIR="${REFERENCE_AUDIT_CACHE:-$PWD/reference_database_audit/uhgg_v2.0.2}"
BASE_URL="${UHGG_AUDIT_BASE_URL:-https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2}"

mkdir -p "$CACHE_DIR"
CACHE_DIR="$(cd "$CACHE_DIR" && pwd -P)"

download() {
  local name="$1"
  local destination="$CACHE_DIR/$name"
  local url="$BASE_URL/$name"

  if [[ -s "$destination" ]]; then
    echo "[OK] Reusing: $destination"
    return
  fi

  echo "[INFO] Downloading: $url"
  echo "[INFO] Destination: $destination"
  if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 5 --retry-delay 5 \
      --continue-at - --output "$destination.part" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget --continue --tries=5 --waitretry=5 \
      --output-document="$destination.part" "$url"
  else
    echo "[ERROR] Neither curl nor wget is available." >&2
    exit 1
  fi
  [[ -s "$destination.part" ]] || {
    echo "[ERROR] Download is empty: $destination.part" >&2
    exit 1
  }
  mv "$destination.part" "$destination"
  echo "[OK] Downloaded: $destination"
}

cat <<EOF
[INFO] Preparing UHGG v2.0.2 reference-audit assets.
[INFO] Expected download is approximately 2.3 GB:
       - genomes-all_metadata.tsv (approximately 110 MB)
       - all_genomes.msh (approximately 2.2 GB)
[INFO] Interrupted downloads are retained as .part files and resumed.
EOF

download genomes-all_metadata.tsv
download all_genomes.msh

header="$(head -n 1 "$CACHE_DIR/genomes-all_metadata.tsv")"
for column in Genome Genome_accession Species_rep Lineage; do
  if [[ $'\t'"$header"$'\t' != *$'\t'"$column"$'\t'* ]]; then
    echo "[ERROR] UHGG metadata lacks required column: $column" >&2
    exit 1
  fi
done

{
  echo "release=UHGG v2.0.2"
  echo "source=$BASE_URL"
  echo "prepared_at=$(date --iso-8601=seconds)"
  echo "cache_dir=$CACHE_DIR"
  sha256sum "$CACHE_DIR/genomes-all_metadata.tsv"
  # Hashing the 2.2-GB sketch is intentional and occurs only after download.
  sha256sum "$CACHE_DIR/all_genomes.msh"
} > "$CACHE_DIR/ASSET_PROVENANCE.txt"

echo "[PASS] Reference-audit assets are ready."
echo "UHGG_ROOT=$CACHE_DIR"
