#!/usr/bin/env bash
# Copy immutable finalized spike pools to fast site-local/shared scratch.
set -euo pipefail
IFS=$'\n\t'

usage() {
  echo "Usage: $0 --source POOLS_DIR --destination SSD_POOLS_DIR" >&2
}
source_dir=""; destination=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source) source_dir="$2"; shift 2 ;;
    --destination) destination="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
[[ -d "$source_dir" && -s "$source_dir/pool_files.sha256" ]] || {
  echo "[ERROR] Source is not a finalized pool directory: $source_dir" >&2
  exit 1
}
[[ -n "$destination" && "$destination" == /* ]] || {
  echo "[ERROR] Destination must be an absolute path" >&2
  exit 2
}
source_dir="$(cd "$source_dir" && pwd -P)"
mkdir -p "$destination"
destination="$(cd "$destination" && pwd -P)"
[[ "$source_dir" != "$destination" ]] || { echo "[ERROR] Source and destination are identical" >&2; exit 1; }

mapfile -t pool_files < <(awk '{print $2}' "$source_dir/pool_files.sha256" | sed 's|^\./||')
(( ${#pool_files[@]} > 0 )) || { echo "[ERROR] Empty checksum manifest" >&2; exit 1; }
for file in "${pool_files[@]}"; do
  [[ "$file" != */* && -s "$source_dir/$file" ]] || {
    echo "[ERROR] Unsafe or missing pool file from manifest: $file" >&2
    exit 1
  }
done

echo "[INFO] Staging ${#pool_files[@]} immutable pool files"
rsync -a --partial --info=progress2 \
  "${pool_files[@]/#/$source_dir/}" \
  "$source_dir/pool_files.sha256" \
  "$source_dir/pool_pair_counts.tsv" \
  "$source_dir/pool_pair_counts.tsv.sha256" \
  "$destination/"

(
  cd "$destination"
  sha256sum -c pool_files.sha256
)
expected_index_sha="$(awk 'NR == 1 {print $1}' "$source_dir/pool_pair_counts.tsv.sha256")"
actual_index_sha="$(sha256sum "$destination/pool_pair_counts.tsv" | awk '{print $1}')"
[[ "$expected_index_sha" == "$actual_index_sha" ]] || {
  echo "[ERROR] Staged pool pair-count index checksum mismatch" >&2
  exit 1
}
(
  cd "$destination"
  sha256sum pool_pair_counts.tsv > pool_pair_counts.tsv.sha256
)
printf 'source\t%s\nstaged_utc\t%s\n' "$source_dir" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  > "$destination/STAGED_SUCCESS"
echo "[PASS] Verified SSD pool cache: $destination"
