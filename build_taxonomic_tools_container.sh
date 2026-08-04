#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT"
DEF="${DEF:-installation/Apptainer.def}"
: "${SIF:?Set SIF to the desired final taxonomic-tools image path}"
FINAL_SIF="$SIF"
BUILD_TMPDIR="${BUILD_TMPDIR:-${TMPDIR:-/tmp}}"
mkdir -p "$BUILD_TMPDIR" "$(dirname "$FINAL_SIF")"
TMP_SIF="$(mktemp "$BUILD_TMPDIR/taxonomic_tools_build.XXXXXX.sif")"
rm -f "$TMP_SIF"

cleanup() {
  [[ ! -e "$TMP_SIF" ]] || rm -f "$TMP_SIF"
}
trap cleanup EXIT

[[ -s "$DEF" ]] || { echo "[ERROR] Missing definition: $DEF" >&2; exit 1; }
if [[ -s "$FINAL_SIF" && "${ALLOW_OVERWRITE:-false}" != "true" ]]; then
  echo "[ERROR] Final image exists: $FINAL_SIF" >&2
  echo "        Use a new versioned filename, or set ALLOW_OVERWRITE=true intentionally." >&2
  exit 1
fi

echo "[INFO] Building upstream image under local temporary storage: $TMP_SIF"
if [[ "${BUILD_WITHOUT_FAKEROOT:-false}" == "true" ]]; then
  echo "[INFO] Using the site-default builder (BUILD_WITHOUT_FAKEROOT=true)"
  apptainer build "$TMP_SIF" "$DEF"
else
  apptainer build --fakeroot "$TMP_SIF" "$DEF"
fi
[[ -s "$TMP_SIF" ]] || { echo "[ERROR] Build produced no image: $TMP_SIF" >&2; exit 1; }

# Some clusters allow building in /tmp but execute images only from an approved
# directory. Install first, then run every verification from the final path.
install -m 0644 "$TMP_SIF" "$FINAL_SIF"
rm -f "$TMP_SIF"

echo "[INFO] Verifying installed upstream image: $FINAL_SIF"
apptainer exec --cleanenv "$FINAL_SIF" /opt/verify_taxonomic_tools.sh
sha256sum "$FINAL_SIF" > "${FINAL_SIF}.sha256"
apptainer inspect "$FINAL_SIF" > "${FINAL_SIF}.inspect.txt"
ls -lh "$FINAL_SIF"
echo "[PASS] Built, installed, and verified: $FINAL_SIF"
