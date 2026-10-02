# Combined pre-production checks

No automatic production authorization. One submission snapshots committed
code, prepares plans, runs two independent arrays, then collects a combined
report after both arrays terminate. Unrelated live worktree edits are excluded.

## Work delivered

- Compact DA3 draft plan: 100 allocations, full seven-dose uniform and
  five-dose partial grids, variable/partial-variable/null arms, four n values,
  both backgrounds/tools and three cohorts: 120,000 contexts in 4,800 batches.
  Only 48 representative batches (1,200 contexts) run in this check.
- Compact workers retain compressed target-level p/adjusted values and native
  source hashes, all-family discovery counts and timing. Full non-target p
  tables are reproducible from plan/input/source but not saved by default.
  Metadata/exposure/dose assignments remain in the immutable plan. This is not
  a full-feature archival package or final production retry/resume runner.
- DA1: twelve actual pilot native matrices/designs checked against pinned
  MaAsLin2 age/sex-adjusted backend, plus 36 synthetic conditional-design null
  stress tasks (Gaussian, skewed and group-heteroskedastic).
- DA3: 24 partial-null synthetic tasks, n=5/10/15/20, family sizes 474/3471,
  Gaussian/correlated/skewed outcomes. Ten features have group shift +1;
  others have unchanged distributions. 200 repetitions/task, exact at n=5,
  otherwise 9,999 MC draws. False-discovery counts refer only to known null
  features. This is a stress check, not a proof of arbitrary strong control.

Both arrays default to 15 concurrent tasks each (30 combined). Preparation
records image SHA and source commit. A failed preparation allows tasks to fail
and the after-any collector to report incompleteness, rather than leaving an
unsatisfiable after-ok dependency indefinitely. No failed task is treated as
a nondiscovery. Report status always requires scientific/resource review.

## Remaining boundaries

Actual MaAsLin2 is not installed locally; backend validation deliberately
requires cluster image version 1.18.0. Local tests exercise real native parsing,
exact inference, compact retention, partial-null simulation, covariate-design
synthetic nulls and failure reporting. Review observed error rates with their
Monte Carlo intervals, not a software PASS alone.

The 100-allocation draft is NOT submitted for production. Shard scheduling
limits, achieved-dose ledger linkage, explicit retry policy, non-target detail
retention and reconciled final protocol still need review. Simulated effects
are on transformed scale, not artificial observed disease profiles. Synthetic
covariate error rates do not prove causal identification in observational DA1.

Compact canary timing includes R startup and shared-file reads; extrapolate
with uncertainty and include scheduler/storage overhead. No significance-based
choice of repeat budget or method is allowed.
