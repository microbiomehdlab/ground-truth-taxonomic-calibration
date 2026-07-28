#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage(){
  cat <<'EOF'
Usage:
  build_samples_spiked_all.sh \
    --samples samples.tsv \
    --independent-root /path/to/WORK/independent \
    --out samples_spiked_all.tsv \
    [--include-original true|false]

What it does:
  - Writes the header from the original samples.tsv
  - Optionally appends original (unspiked) samples
  - Appends ALL spiked entries found under:
      <independent-root>/<LABEL>/logs/spiked_manifest_parts/*.tsv

Assumptions:
  - samples.tsv has header and at least: sample_id, fastq1, fastq2
  - spiked_manifest_parts rows already preserve extra columns (e.g., original_id, Target_Condition)

EOF
}

SAMPLES=""
INDEP=""
OUT=""
INCLUDE_ORIG="true"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2;;
    --independent-root) INDEP="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --include-original) INCLUDE_ORIG="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$SAMPLES" && -n "$INDEP" && -n "$OUT" ]] || { usage; exit 2; }
[[ -f "$SAMPLES" ]] || { echo "[ERROR] Missing samples file: $SAMPLES" >&2; exit 1; }
[[ -d "$INDEP" ]] || { echo "[ERROR] Missing independent root dir: $INDEP" >&2; exit 1; }

case "$INCLUDE_ORIG" in true|false) ;; *) echo "[ERROR] --include-original must be true|false" >&2; exit 2;; esac

tmp="${OUT}.tmp"
: > "$tmp"

# header (exactly as original samples.tsv)
head -n 1 "$SAMPLES" >> "$tmp"

# optionally include original rows
if [[ "$INCLUDE_ORIG" == "true" ]]; then
  tail -n +2 "$SAMPLES" >> "$tmp"
fi

# append all spiked part files, deterministic order
found=0
for label_dir in "$INDEP"/*; do
  [[ -d "$label_dir" ]] || continue
  parts="$label_dir/logs/spiked_manifest_parts"
  [[ -d "$parts" ]] || continue

  # append in sorted order
  while IFS= read -r f; do
    [[ -s "$f" ]] || continue
    cat "$f" >> "$tmp"
    found=1
  done < <(find "$parts" -type f -name "*.tsv" | sort)
done

[[ "$found" -eq 1 ]] || echo "[WARN] No spiked manifest parts found under $INDEP" >&2

mv -f "$tmp" "$OUT"
echo "[OK] Wrote $OUT"