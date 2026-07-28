#!/usr/bin/env bash
set -euo pipefail

DESTINATION="${DESTINATION:-$PWD}"
SOURCE_DATA_ROOT="${SOURCE_DATA_ROOT:-}"
ALLOW_OVERWRITE="${ALLOW_OVERWRITE:-false}"

if [[ -z "$SOURCE_DATA_ROOT" ]]; then
  echo "[ERROR] Set SOURCE_DATA_ROOT to the directory containing the required inputs." >&2
  exit 1
fi

SOURCE_DATA_ROOT="$(cd "$SOURCE_DATA_ROOT" && pwd -P)"
DESTINATION="$(cd "$DESTINATION" && pwd -P)"

if [[ "$SOURCE_DATA_ROOT" == "$DESTINATION" ]]; then
  echo "[ERROR] SOURCE_DATA_ROOT and DESTINATION are the same directory." >&2
  exit 1
fi

input_files=(
  samples_spiked_all_independent.tsv
  samples_spiked_all_community.tsv
  spike_panel.tsv
  spike_taxon_aliases.csv
  metadata_w_study.tsv
  metaphlan4_merged_unspecified.csv
  kraken2_bracken_merged_unspecified.csv
)

for file in "${input_files[@]}"; do
  [[ -s "$SOURCE_DATA_ROOT/$file" ]] || {
    echo "[ERROR] Missing or empty source input: $SOURCE_DATA_ROOT/$file" >&2
    exit 1
  }
done

[[ -d "$SOURCE_DATA_ROOT/profile_results" ]] || {
  echo "[ERROR] Missing source directory: $SOURCE_DATA_ROOT/profile_results" >&2
  exit 1
}
find "$SOURCE_DATA_ROOT/profile_results" -type f -print -quit | grep -q . || {
  echo "[ERROR] Source profile_results directory is empty." >&2
  exit 1
}
for file in "${input_files[@]}"; do
  target="$DESTINATION/$file"
  if [[ -e "$target" && "$ALLOW_OVERWRITE" != "true" ]]; then
    echo "[ERROR] Destination already exists: $target" >&2
    echo "        Set ALLOW_OVERWRITE=true only to replace staged inputs." >&2
    exit 1
  fi
done
if [[ -e "$DESTINATION/profile_results" && ! -d "$DESTINATION/profile_results" ]]; then
  echo "[ERROR] Destination exists but is not a directory: $DESTINATION/profile_results" >&2
  exit 1
fi

rsync_options=(-a --info=progress2)
if [[ "$ALLOW_OVERWRITE" != "true" ]]; then
  rsync_options+=(--ignore-existing)
fi

for file in "${input_files[@]}"; do
  rsync "${rsync_options[@]}" "$SOURCE_DATA_ROOT/$file" "$DESTINATION/$file"
done
rsync "${rsync_options[@]}" \
  "$SOURCE_DATA_ROOT/profile_results/" \
  "$DESTINATION/profile_results/"
{
  for file in "${input_files[@]}"; do
    sha256sum "$DESTINATION/$file"
  done
  find "$DESTINATION/profile_results" -type f -print0 |
    sort -z |
    xargs -0 sha256sum
} > "$DESTINATION/LOCAL_INPUT_SHA256.txt"

echo "[PASS] Required inputs staged into: $DESTINATION"
du -sh "${input_files[@]/#/$DESTINATION/}" "$DESTINATION/profile_results"
echo "[INFO] Checksums: $DESTINATION/LOCAL_INPUT_SHA256.txt"
