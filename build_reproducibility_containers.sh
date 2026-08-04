#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$ROOT"
: "${UPSTREAM_SIF:?Set the final taxonomic-tools SIF path}"
: "${ANALYSIS_SIF:?Set the final MaAsLin2/Mash SIF path}"

if [[ "${ALLOW_DIRTY_SOURCE:-false}" != "true" ]] &&
   [[ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no 2>/dev/null || true)" ]]; then
  echo "[ERROR] Tracked source files are modified. Commit or restore them before building." >&2
  echo "        Set ALLOW_DIRTY_SOURCE=true only for an explicitly non-publication test." >&2
  exit 1
fi

echo "[STEP 1/2] Build preprocessing/profiling image"
SIF="$UPSTREAM_SIF" \
DEF="${UPSTREAM_DEF:-installation/Apptainer.def}" \
BUILD_TMPDIR="${BUILD_TMPDIR:-${TMPDIR:-/tmp}}" \
ALLOW_OVERWRITE="${ALLOW_OVERWRITE:-false}" \
bash "$ROOT/build_taxonomic_tools_container.sh"

echo "[STEP 2/2] Build downstream analysis image"
SIF="$ANALYSIS_SIF" \
DEF="${ANALYSIS_DEF:-containers/Apptainer.def}" \
BUILD_TMPDIR="${BUILD_TMPDIR:-${TMPDIR:-/tmp}}" \
ALLOW_OVERWRITE="${ALLOW_OVERWRITE:-false}" \
bash "$ROOT/build_crc_spike_maaslin2_container.sh"

manifest="${CONTAINER_MANIFEST:-$ROOT/work/container_build/container_build_manifest.tsv}"
mkdir -p "$(dirname "$manifest")"
hash_inputs() {
  sha256sum "$@" | sha256sum | awk '{print $1}'
}
git_commit="$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || printf 'not-a-git-checkout')"
{
  printf 'role\tpath\tsha256\tdefinition\tdefinition_sha256\trecipe_inputs_sha256\tgit_commit\n'
  printf 'upstream\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$UPSTREAM_SIF" "$(sha256sum "$UPSTREAM_SIF" | awk '{print $1}')" \
    "${UPSTREAM_DEF:-installation/Apptainer.def}" \
    "$(sha256sum "${UPSTREAM_DEF:-installation/Apptainer.def}" | awk '{print $1}')" \
    "$(hash_inputs "${UPSTREAM_DEF:-installation/Apptainer.def}" installation/verify_taxonomic_tools.sh)" \
    "$git_commit"
  printf 'analysis\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$ANALYSIS_SIF" "$(sha256sum "$ANALYSIS_SIF" | awk '{print $1}')" \
    "${ANALYSIS_DEF:-containers/Apptainer.def}" \
    "$(sha256sum "${ANALYSIS_DEF:-containers/Apptainer.def}" | awk '{print $1}')" \
    "$(hash_inputs "${ANALYSIS_DEF:-containers/Apptainer.def}" containers/environment.yml containers/verify_crc_spike_environment.R)" \
    "$git_commit"
} > "$manifest"
apptainer version > "${manifest%.tsv}.apptainer_version.txt"

echo "[PASS] Both source-defined images were rebuilt and verified"
echo "[INFO] Build manifest: $manifest"
