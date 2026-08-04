#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT"
DEF="${DEF:-containers/Apptainer.def}"
: "${SIF:?Set SIF to the desired final analysis-container path}"
FINAL_SIF="$SIF"
BUILD_TMPDIR="${BUILD_TMPDIR:-${TMPDIR:-/tmp}}"
mkdir -p "$BUILD_TMPDIR"
TMP_SIF="$(mktemp "$BUILD_TMPDIR/crc_spike_maaslin2_build.XXXXXX.sif")"
rm -f "$TMP_SIF"

cleanup() {
  if [[ -e "$TMP_SIF" ]]; then
    rm -f "$TMP_SIF"
  fi
}
trap cleanup EXIT

mkdir -p "$(dirname "$FINAL_SIF")"

echo "[INFO] Definition: $DEF"
echo "[INFO] Temporary build SIF: $TMP_SIF"
echo "[INFO] Final SIF: $FINAL_SIF"

if [[ -s "$FINAL_SIF" && "${ALLOW_OVERWRITE:-false}" != "true" ]]; then
  echo "[ERROR] Final SIF already exists and is non-empty: $FINAL_SIF" >&2
  echo "        Set ALLOW_OVERWRITE=true only if you intentionally want to replace it." >&2
  exit 1
fi

echo "[INFO] Building temporary image..."
if [[ "${BUILD_WITHOUT_FAKEROOT:-false}" == "true" ]]; then
  echo "[INFO] Using the site-default builder (BUILD_WITHOUT_FAKEROOT=true)"
  apptainer build "$TMP_SIF" "$DEF"
else
  apptainer build --fakeroot "$TMP_SIF" "$DEF"
fi

echo "[INFO] Checking temporary image..."
if [[ ! -s "$TMP_SIF" ]]; then
  echo "[ERROR] Temporary SIF is missing or empty: $TMP_SIF" >&2
  exit 1
fi

ls -lh "$TMP_SIF"
file "$TMP_SIF"

echo "[INFO] Installing SIF into the Apptainer-approved image directory..."
# Build under local scratch because fakeroot may be unable to create its output
# directly on NFS/BeeGFS even when the calling user can write there.
install -m 0644 "$TMP_SIF" "$FINAL_SIF"
rm -f "$TMP_SIF"

echo "[INFO] Verifying R environment from the installed image..."
# Some clusters permit building under /tmp but restrict image execution to an
# administrator-approved directory such as /mnt/beegfs/apptainer/images.
apptainer exec --cleanenv "$FINAL_SIF" \
  Rscript /opt/verify_crc_spike_environment.R
apptainer exec --cleanenv "$FINAL_SIF" mash --version
sha256sum "$FINAL_SIF" > "${FINAL_SIF}.sha256"
apptainer inspect "$FINAL_SIF" > "${FINAL_SIF}.inspect.txt"

ls -lh "$FINAL_SIF"
file "$FINAL_SIF"

echo "[PASS] Built and installed: $FINAL_SIF"
