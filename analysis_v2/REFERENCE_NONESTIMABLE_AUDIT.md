# Saved-input audit for reference pilot NE cases

Purpose: distinguish all-zero/nonzero constant MetaPhlAn inputs from perfect
fixed-group Bracken fits using the checksum-verified cache actually supplied to
the worker. No live FASTQs/profiles, new profiling, model refitting, p-values or
production authorization. The private downloaded report lacks these caches, so
execution must occur on the cluster; local fixtures are not real-data findings.

```bash
export PROJECT="$PWD"
export REFERENCE_ROOT="$PWD/work/reference_response_pilot_20261005T202107Z"
export REFERENCE_NE_ROOT="$PWD/work/reference_ne_audit_$(date -u +%Y%m%dT%H%M%SZ)"
REFERENCE_NE_JOB="$(sbatch --parsable --export=ALL analysis_v2/audit_reference_nonestimable.sbatch)"
printf 'Job: %s\nOutput: %s\n' "$REFERENCE_NE_JOB" "$REFERENCE_NE_ROOT"
```

Read `status.json` after job completion and verify `SHA256SUMS`. Download this
small audit bundle, not the full analysis cache. Expected real-data scope is
439 observed NE target/context rows (148 Bracken, 291 MetaPhlAn), with repeated
people across contexts—not 439 independent biological failures.

`input_summary.tsv` preserves saved model status, raw native range/zero counts,
group means, transformed range, arithmetic residual check and numerical tolerance.
`input_vectors.tsv` contains every person's group, tested dose, raw abundance,
transformed value and source path for these targets. `provenance.json` binds the
audit to the plan/report manifests and audit script. Neither input bundle nor
completed model results are changed. A new sibling output is required.

The diagnostic uses the worker's log2(1+a/1e-8) transform and tolerance
100 × machine epsilon × max(1,max(abs(y))). Constant range is checked first;
otherwise residuals from the two group means detect perfect fixed-group fits.
This arithmetic check is not an independent MaAsLin fit; disagreement with saved
status stops the audit for review. Absence from a saved profile maps to zero,
exactly as in the worker. All-zero profile input establishes analytical zero
reporting, not absence of organisms or absence of classified reads.

Tests: all-zero/nonzero constants, perfect separation, variable inputs, malformed
vectors, complete sealed-cache audit, refusal to overwrite and checksum tampering.
