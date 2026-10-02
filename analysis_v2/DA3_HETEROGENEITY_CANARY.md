# DA3 heterogeneity engineering canary

This validates assembly and measures inference costs; it does not authorize
production, validate partial-null error control or select a method by its hits.
The proposal is in `PRODUCTION_EXPERIMENT_MATRIX_20261002.md`.

384 contexts: three cohorts, adenoma/CRC, both tools, n=5/10/15/20,
and eight scenarios: uniform at two anchors, partial at 25/50/75%, variable,
partial-variable at 50%, and no added signal. One allocation per stratum.
Full production dose/repetition coverage is deliberately not generated.

Exposed cases use actual existing profiles; unexposed cases and controls use
originals. Dose assignment and exposure ranks are shared across tools and
independent of measured abundance. Half-up exposure rounding is explicit.
Variable amounts use the fixed five-dose grid. Metadata record nominal amounts;
achieved amounts still need the canonical/design ledger for production.

The planner verifies the existing inventory, freezes its/code hashes and
selects the baseline frozen family. Workers verify plan/source profiles and
code hashes, retain native matrices and run the existing exact/MC candidate
plus parametric comparison. Exact n=5; MC 9,999 permutations otherwise.
One fixed MC seed is an engineering choice, not a production seed policy.

Fresh output directories only; failed inference has no SUCCESS. Collector
reports every missing/failed context and never authorizes production.
Timing is separate assembly and inference cost, not yet a wall-time forecast
including scheduling, storage and parallel contention. All full-family feature
results are retained for this small canary, not proposed for millions of contexts.

Do not modify/pull the canary source files while workers run: hash checks
reject drift. A production runner must use immutable source snapshots.
The canary is bounded and independently parallelizable but intentionally
does not claim production retry/resume semantics.

Local verification: three fixture tests cover full 384-context planning,
deterministic assignments/rounding, no biological overlap, actual native
profile parsing and exact inference, tamper rejection and failure-aware
collection. Large-n MC and exact functions have their existing separate tests.
