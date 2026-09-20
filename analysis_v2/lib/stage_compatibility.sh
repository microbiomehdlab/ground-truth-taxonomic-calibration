#!/usr/bin/env bash
# Fail-closed reuse gate for completed pipeline stage directories.
#
# A `SUCCESS` marker only proves the stage finished under the code of the day.
# When a producer gains columns, an older completed directory is still marked
# SUCCESS but can no longer feed the current consumers, so reusing it fails the
# run downstream instead of at the reuse decision. This library asks
# `scripts/check_stage_schema.py` whether a directory matches the current
# schema, and quarantines it under `failed_attempts/` when it does not.
#
# Header validation is exact tab-separated parsing inside that helper; never an
# unquoted grep, which would match `primary_reference_type` when
# `primary_row_reference_type` is required.
#
# Sourced by analysis_v2/run_geff_propagation_development.sh.

STAGE_COMPATIBILITY_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
STAGE_SCHEMA_CHECKER="$STAGE_COMPATIBILITY_LIB_DIR/../scripts/check_stage_schema.py"

# Use the container interpreter when the caller provides one, so the runner and
# a bare test harness both work without extra plumbing.
stage_compatibility_python() {
  if declare -F analysis_python >/dev/null 2>&1; then
    analysis_python "$@"
  else
    python3 "$@"
  fi
}

# stage_schema_reason STAGE DIR
# Prints the first incompatibility reason and returns 1; returns 0 when reusable.
stage_schema_reason() {
  local stage="$1" dir="$2"
  stage_compatibility_python "$STAGE_SCHEMA_CHECKER" \
    --stage "$stage" --directory "$dir" --quiet
}

# reuse_or_quarantine STAGE DIR RUN_ROOT STAMP LABEL
#
# Returns 0 when DIR is complete, schema-compatible and may be reused.
# Otherwise DIR (if present) is moved whole to
# RUN_ROOT/failed_attempts/<basename>_<STAMP> and 1 is returned, so the caller
# regenerates it. Inputs of the stage are never touched.
reuse_or_quarantine() {
  local stage="$1" dir="$2" run_root="$3" stamp="$4" label="$5"
  if [[ ! -e "$dir" ]]; then
    return 1
  fi
  local reason
  if reason="$(stage_schema_reason "$stage" "$dir")"; then
    echo "[REUSE] $label: $dir"
    return 0
  fi
  mkdir -p "$run_root/failed_attempts"
  local quarantine="$run_root/failed_attempts/$(basename "$dir")_${stamp}"
  mv "$dir" "$quarantine"
  echo "[QUARANTINE] stale or incompatible $label moved under failed_attempts"
  echo "[QUARANTINE]   reason: ${reason:-unknown}"
  echo "[QUARANTINE]   kept at: $quarantine"
  return 1
}
