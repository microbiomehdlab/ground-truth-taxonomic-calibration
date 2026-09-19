#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT"

DEF="${DEF:-containers/Apptainer.analysis_v2.def}"
: "${SIF:?Set SIF to a new, versioned final image path}"
FINAL_SIF="$SIF"
BUILD_TMPDIR="${BUILD_TMPDIR:-${TMPDIR:-/tmp}}"

[[ "$FINAL_SIF" == *.sif ]] || { echo "[ERROR] SIF must end in .sif" >&2; exit 1; }
[[ ! -e "$FINAL_SIF" ]] || {
  echo "[ERROR] Refusing to overwrite existing image: $FINAL_SIF" >&2
  exit 1
}
[[ -s "$DEF" ]] || { echo "[ERROR] Missing definition: $DEF" >&2; exit 1; }

mkdir -p "$BUILD_TMPDIR" "$(dirname "$FINAL_SIF")"
TMP_SIF="$(mktemp "$BUILD_TMPDIR/ground_truth_analysis_v2.XXXXXX.sif")"
rm -f "$TMP_SIF"

cleanup() {
  [[ ! -e "$TMP_SIF" ]] || rm -f "$TMP_SIF"
}
trap cleanup EXIT

echo "[INFO] Definition: $DEF"
echo "[INFO] Temporary image: $TMP_SIF"
echo "[INFO] Final image: $FINAL_SIF"

if [[ "${BUILD_WITHOUT_FAKEROOT:-false}" == "true" ]]; then
  apptainer build "$TMP_SIF" "$DEF"
else
  apptainer build --fakeroot "$TMP_SIF" "$DEF"
fi

[[ -s "$TMP_SIF" ]] || { echo "[ERROR] Built image is empty" >&2; exit 1; }
install -m 0644 "$TMP_SIF" "$FINAL_SIF"
rm -f "$TMP_SIF"

apptainer test "$FINAL_SIF"
apptainer exec --cleanenv "$FINAL_SIF" \
  conda list -n ground-truth-analysis-v2 --explicit \
  > "${FINAL_SIF}.conda-explicit.txt"
apptainer exec --cleanenv "$FINAL_SIF" \
  Rscript -e 'sessioninfo::session_info()' \
  > "${FINAL_SIF}.R-session-info.txt"
apptainer inspect "$FINAL_SIF" > "${FINAL_SIF}.inspect.txt"
sha256sum "$FINAL_SIF" > "${FINAL_SIF}.sha256"
sha256sum \
  "$DEF" \
  containers/environment_analysis_v2.yml \
  containers/verify_ground_truth_analysis_v2.py \
  containers/verify_ground_truth_analysis_v2.R \
  > "${FINAL_SIF}.source-files.sha256"

echo "[PASS] Built and verified: $FINAL_SIF"
cat "${FINAL_SIF}.sha256"
