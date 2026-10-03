#!/usr/bin/env bash
# Read-only audit/export of existing clinical fits. Never submit new model fits.
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA_PILOT_PLAN:?}" "${DA_PILOT_RESULTS:?}" "${DA1_EXPORT_ROOT:?}"
cd "$PROJECT"
test ! -e "$DA1_EXPORT_ROOT"
test -s "$ANALYSIS_SIF"
test -s "$DA_PILOT_PLAN/analysis_image.sha256"
test -s "$DA_PILOT_PLAN/contexts.tsv"
test -d "$DA_PILOT_RESULTS"
mkdir -p "$DA1_EXPORT_ROOT/source"
git archive HEAD | tar -x -C "$DA1_EXPORT_ROOT/source"
git rev-parse HEAD > "$DA1_EXPORT_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$DA1_EXPORT_ROOT/analysis_image.sha256"
(cd "$DA1_EXPORT_ROOT/source" && find . -type f ! -name SOURCE_SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SOURCE_SHA256SUMS)
export PROJECT ANALYSIS_SIF DA_PILOT_PLAN DA_PILOT_RESULTS DA1_EXPORT_ROOT
job="$(sbatch --parsable --export=ALL "$DA1_EXPORT_ROOT/source/analysis_v2/run_da1_clinical_export.sbatch")"
printf 'job\t%s\n' "$job" | tee "$DA1_EXPORT_ROOT/jobs.tsv"
printf 'root\t%s\n' "$DA1_EXPORT_ROOT"
