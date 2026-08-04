#!/usr/bin/env bash
# One-sample preprocessing validation runner. It never deletes scratch inputs.
set -euo pipefail
IFS=$'\n\t'

: "${YACHIDA_ENV:?set YACHIDA_ENV to the local ignored configuration}"
: "${RECEIPT:?provided by stream_sample.py}"
: "${SAMPLE_WORK:?provided by stream_sample.py}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
bash "$ROOT/datasets/yachida/run_metashotgunprep.sh"

# shellcheck source=/dev/null
source "$SAMPLE_WORK/metashotgunprep_outputs.env"
for file in "$CLEAN_R1" "$CLEAN_R2"; do
  [[ -s "$file" ]] || { echo "[ERROR] Smoke-test cleaned mate missing: $file" >&2; exit 1; }
  gzip -t -- "$file"
done

python3 - "$PERSISTENT_SAMPLE_QC" "$RECEIPT" <<'PY'
import hashlib
import pathlib
import sys

qc_root = pathlib.Path(sys.argv[1]).resolve()
receipt = pathlib.Path(sys.argv[2])
files = sorted(path for path in qc_root.rglob("*") if path.is_file() and path.stat().st_size)
if not files:
    raise SystemExit(f"[ERROR] No persistent QC files found under {qc_root}")
receipt.parent.mkdir(parents=True, exist_ok=True)
with receipt.open("w", encoding="utf-8") as handle:
    handle.write("path\tsha256\tbytes\n")
    for path in files:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            while block := source.read(4 * 1024 * 1024):
                digest.update(block)
        handle.write(f"{path}\t{digest.hexdigest()}\t{path.stat().st_size}\n")
print(f"[OK] Recorded {len(files)} persistent preprocessing outputs: {receipt}")
PY

echo "[PASS] One-sample preprocessing smoke test completed: $SAMPLE_ID"
echo "[KEEP] Cleaned reads retained for inspection under: $SAMPLE_WORK"
