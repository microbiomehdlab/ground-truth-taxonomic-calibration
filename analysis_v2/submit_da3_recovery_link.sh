#!/usr/bin/env bash
set -euo pipefail
: "${PROJECT:?}" "${ANALYSIS_SIF:?}" "${DA3_BATCH_PLAN:?}"
: "${BRACKEN_RECOVERY_ROOT:?}" "${METAPHLAN_RECOVERY_ROOT:?}" "${RECOVERY_LINK_ROOT:?Fresh output root required}"
cd "$PROJECT"
if [[ -z "${DA3_RUN_ROOT:-}" ]]; then
    mapfile -t da3_candidates < <(find "$PROJECT/work" -maxdepth 1 -type d -name 'da3_exploratory_100_*' | sort)
    if [[ "${#da3_candidates[@]}" -ne 1 ]]; then
        printf 'Set DA3_RUN_ROOT explicitly; expected one exploratory run, found %s.\n' "${#da3_candidates[@]}" >&2
        printf '%s\n' "${da3_candidates[@]}" >&2
        exit 1
    fi
    export DA3_RUN_ROOT="${da3_candidates[0]}"
fi
test ! -e "$RECOVERY_LINK_ROOT"
test -s "$ANALYSIS_SIF"
test -s "$DA3_BATCH_PLAN/SHA256SUMS"
test -s "$BRACKEN_RECOVERY_ROOT/feng/bracken_recovery_comparison.tsv"
test -s "$METAPHLAN_RECOVERY_ROOT/feng/endpoints/paired_endpoints.tsv"
python3 - "$DA3_RUN_ROOT/REPORT/status.json" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]))
assert s['status']=='COMPLETE_PENDING_SCIENTIFIC_REVIEW' and not s['failures'], 'DA3 collection incomplete'
assert s['completed_contexts']==120000 and s['completed_batches']==s['expected_batches']==4800, 'Unexpected DA3 coverage'
PY
mkdir -p "$RECOVERY_LINK_ROOT/source"
git archive HEAD | tar -x -C "$RECOVERY_LINK_ROOT/source"
git rev-parse HEAD > "$RECOVERY_LINK_ROOT/source_commit.txt"
sha256sum "$ANALYSIS_SIF" > "$RECOVERY_LINK_ROOT/image.sha256"
export RECOVERY_LINK_SOURCE="$RECOVERY_LINK_ROOT/source"
export RECOVERY_LINK_ROOT ANALYSIS_SIF DA3_BATCH_PLAN DA3_RUN_ROOT BRACKEN_RECOVERY_ROOT METAPHLAN_RECOVERY_ROOT
job="$(sbatch --parsable --export=ALL "$RECOVERY_LINK_SOURCE/analysis_v2/run_da3_recovery_link.sbatch")"
printf 'job\t%s\n' "$job" | tee "$RECOVERY_LINK_ROOT/jobs.tsv"
printf 'root\t%s\n' "$RECOVERY_LINK_ROOT"
