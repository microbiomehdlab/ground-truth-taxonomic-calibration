# Fast parallel candidate comparison — 2026-10-02

User requested DA tests first while biological-realism development proceeds.
No definitive method amendment or production approval is implied.

52 independent array tasks: 12 null strata (all three cohorts, two profilers,
Adenoma and CRC), each with the frozen 1,000 full-pool n=5 allocations and all
252 independent-subset n=5 allocations; and all 40 spiked DA3 n=5 pilot contexts.
Use only saved, checksum-verified native matrices, metadata, families and
allocations. No re-profiling, dose changes or outcome-selected feature filtering.

Compare existing shifted-log mean contrasts/parametric BH with exact two-sided
permutation p-values across every one of the 252 assignments. Report both full
family BH and separately named single-step family-maximum adjusted p-values.
The latter targets family-wise error under the global exchangeable-label null,
not the same FDR criterion as BH; do not claim arbitrary strong error control
under partial alternatives, or general clinical validity. Preserve negative
effects and all species, including non-targets. Record effect estimates.

The statistic is absolute mean difference scaled by label-invariant pooled
sample SD. Unlike residual-SE inference, perfect group separation is a valid
permutation statistic, not an invented infinite t or tiny p. Globally constant
features retain p=1. Floating ties use 100 machine-epsilon relative guards.
Exact enumeration includes the observed assignment; no Monte Carlo zero p.
Two-sided resolution at balanced 5+5 is at best 2/252. Report power limitations
honestly; do not narrow the family afterwards to obtain discoveries.

This specifically addresses the worst observed DA3 null failure. It does NOT
settle DA1 clinical covariates, larger-n production, or DA2 asymmetric mean-null
inference. Paired sign flips require symmetry and cannot repair the synthetic
skewness failure merely by replacing the t reference distribution. Those method
questions can be developed while this batch runs. Do not redefine the main
experiment matrix or adopt this candidate without comparative review.

Source is snapshotted with git archive at submission, so new development and
git pulls do not invalidate running jobs. Each task hashes its snapshot code
and input matrices; image checks run before/after. The after-any collector
reports incomplete tasks rather than silently dropping failed results.

Primary implementation reference for exact permutation enumeration and tie
handling: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html
No SciPy runtime dependency is required by this R implementation.
