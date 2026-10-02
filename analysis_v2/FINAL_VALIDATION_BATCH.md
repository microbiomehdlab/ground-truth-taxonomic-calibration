# Consolidated validation batch — 2026-10-02

One submission queues 40 independent tasks at up to 30 concurrent tasks, then
an after-any collector. No cluster access by Codex, no primary data mutation,
no production launch. Existing unrelated working-tree edits remain untouched.

1. Twelve paired DA2 pilot contexts: pseudocounts 1e-9, 1e-8 and 1e-7, with
   unchanged native inputs/pairs/families. Retain every result, not only targets.
2. Twelve DA3 strata (three cohorts, two profilers, Adenoma and CRC): 1,000
   baseline-only randomized allocations for each n=5,10,15,20 plus all 252
   matched-subset allocations at n=5. Total 51,024 family-level null fits.
   Reuse frozen inventory allocations/families, never choose by result.
   Fast group-only equal-variance LM computation must agree with the actual
   pinned MaAsLin2 wrapper for the first allocation in EACH pool/n stratum.
   This is an acceleration for validation, not a replacement of DA3's primary
   production software. No silent warnings/failed contexts discarded.
3. Sixteen paired synthetic stress contexts: n=10,42,47,67, 3,471 features,
   1,000 replicates for independent Gaussian, correlated Gaussian, sparse
   symmetric and skewed zero-mean exponential differences. The last scenario
   deliberately tests non-normal small-sample robustness; do not interpret it
   as the actual microbiome null distribution. Independent Monte Carlo draws
   preserve full BH families and explicit numerical degeneracy handling.

Null review flags: any-discovery rate >0.075 or Monte Carlo CI lower bound >0.05.
These are prespecified engineering-review flags, not scientific significance
tests. Under a global null, any two-sided discovery is the relevant error event;
positive counts are additional. Exhaustive subset summaries have no Monte Carlo
CI; full-pool intervals are conditional on the fixed baseline profiles.

## Existing evidence recorded

Sign-flip job 3097937: all 12 contexts completed and checksums verified.
Conditional any-discovery rates 0–0.009. Zeller–MetaPhlAn community 0.009,
mean discoveries 0.196; indicates clustered discoveries in rare sign patterns.
This is reassuring conditionally, not biological FDR certification.
See DA2_PILOT_REVIEW_20261002.md for the preceding target findings.

## Handoff

`bash analysis_v2/submit_final_validation.sh` enforces a real-container backend
fixture before submission. Supply PROJECT, ANALYSIS_SIF, DA2_PAIRED_RESULTS,
DA_INVENTORY_ROOT. Output has jobs.tsv, per-task hashed diagnostics and REPORT.
The collector reports failures even when individual tasks fail. It never
labels production authorized; review the combined summary, pseudocount target
changes and null flags before freezing final policy. Full-dose production
planning/aggregation still needs an explicit final configuration; this batch
does not claim that code is already complete.
