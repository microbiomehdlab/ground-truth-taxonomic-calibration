# Recovery assessment: decisions and implementation checkpoint

## Bracken: agreed primary measurement

Use the all-input, count-based recovery assessment as the primary quantitative
spike-recovery endpoint. Keep the previous native-reported-fraction calculation
as a labelled sensitivity analysis. Do not replace the native abundance tables
used for differential-abundance analysis with these recovery endpoints.

For original input pairs R, total inserted pairs N, target inserted pairs n_t,
baseline Bracken estimated count B_t and spiked estimated count O_t:

- Observed all-input abundance: O_t / (R + N).
- Expected abundance: (B_t + n_t) / (R + N).
- Recovered added signal: (O_t - B_t) / (R + N).
- Recovery ratio: (O_t - B_t) / n_t.

The denominator includes unclassified input pairs. This is an assessment of
changes in estimated counts, not direct tracking of implanted read origins or
absolute microbial cell abundance. Baseline reassignment can influence the
response. Community member counts are reconstructed from the frozen ordered
panel and deterministic integer allocator, not independently read-origin
verified. See BRACKEN_ALL_INPUT_RECOVERY.md for the full contract.

## Real-data evidence reviewed on 2026-10-01

Local immutable report reviewed:
../bracken_report/bracken_recovery_report_20261001T131003Z

Cluster report source commit: 0dd55d4 (report workflow).
Recovery inputs: work/bracken_all_input_recovery_20260930T153115Z.
The report's SHA256SUMS verified locally. Its tables/summary.json records a
non-fixture PASS_SUMMARY_TABLES result, 41,170 target endpoints and 1,170
cohort/population/condition/target/dose summary cells.

All 1,170 all-input cell median recovery ratios are between 0.8 and 1.2
(approximately 0.809–1.051). This describes cell medians, NOT every sample.
Individual endpoint ratios include zero and values above two. Only one target
endpoint has a native/count nondetection; detection alone is not quantitative
accuracy.

Median across context-cell median all-input recovery ratios (descriptive,
not pooled raw-sample estimates or inferential results):

| Target | Ratio |
| --- | ---: |
| Bfrag | 0.992 |
| Csym | 1.031 |
| Dpne | 1.008 |
| Fnuc | 0.907 |
| Hhat | 1.027 |
| Pana | 0.976 |
| Pint | 0.834 |
| Pmic | 0.937 |
| Porp | 0.815 |
| Psto | 0.927 |

Interpretation: recovery is generally stable across doses and clinical
backgrounds, with persistent taxon-specific under-recovery. The new reference
does not improve the profiler itself and does not reduce numerical error for
every target. Native-fraction rounding is particularly relevant at the weakest
community spike; near-perfect native ratios should not be interpreted as proof
of perfect recovery. A causal decomposition of denominator and rounding effects
has not been established by this summary.

Actual retained nominal total doses:
community 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1;
independent 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05.
Community dose is total panel dose, not each member's dose.

Suggested manuscript wording:
“Bracken recovered implanted signals consistently across cohorts and clinical
backgrounds, while all-input read accounting revealed modest, taxon-specific
under-recovery.” Include the sample distributions and reconstruction caveat.
No upstream rerun is indicated by this report.

## MetaPhlAn: existing code, not yet a validated three-cohort result

Already implemented:

- scripts/build_target_genome_sizes.py: implanted assembly lengths.
- scripts/extract_metaphlan_genome_sizes.py and
  scripts/compute_effective_genome_size.py: reference-size inputs.
- scripts/derive_paired_endpoints.py: explicit genome_equivalent primary and
  read_proportional sensitivity modes, using unchanged observed abundances.
- lib/metaphlan_reference.sh: required reference-table checks and separate
  output directories for both modes.
- scripts/compare_target_recovery_references.py: paired reference comparison.
- scripts/audit_metaphlan_genome_size_residual.py: existing Yachida-oriented
  residual diagnostic; not evidence of three-cohort replication.

The primary expectation converts implanted read fractions f_it to
q_it = f_it * G_eff,i / G_t, then uses D = 1 - F + sum_t(q_it).
Expected abundance is [(1 - F) * baseline abundance + q_it] / D.
This changes the expected reference, not the raw MetaPhlAn output. Effective
genome size is an estimated input whose provenance and uncertainty require
review; this is not a guarantee that all residual profiler bias is removed.

Focused verification performed locally on this date:
test_metaphlan_reference_runner.py PASS;
test_target_recovery_reference_comparison.py: 9 tests, OK.
These tests do not establish complete runner readiness or real-data coverage.

Required next work before a full downstream submission:

1. Resolve the confirmed shared-runner directory conflict: the runner creates
   RUN_ROOT/endpoints before derive_endpoints_with_references rejects an
   existing endpoint output. Add an end-to-end regression test.
2. Make preflight validate the required reference tables; its current early
   PREFLIGHT_ONLY exit precedes endpoint reference validation.
3. Verify exact assembly sizes and sample-wide effective-size provenance and
   coverage for all 511 samples, including both analysis populations.
4. Produce primary and sensitivity recovery comparisons for all three cohorts
   and inspect their distributions before making profiler-ranking claims.
5. Update the written DA protocol to the agreed MaAsLin2-primary decision,
   with custom tests as sensitivity, before implementing the revised DA plan.
6. Check that recovery-scale consumers select the intended endpoint columns,
   and that DA still consumes native observed abundance rather than recovery.

No cluster submission, upstream output change, or statistical model run was
performed by this documentation checkpoint.

## 2026-10-01 — Execution-gate fixes implemented locally

The shared runner and both maintained cohort-specific runners no longer
precreate the endpoint directory. All three call require_reference_inputs before
their preflight-only exit. The helper now invokes the read-only validator
validate_metaphlan_reference_inputs.py, reusing endpoint-stage table loaders
and population-specific/sample-wide lookup rules. It rejects nonfinite and
nonpositive sizes, duplicate keys, missing targets and missing population scopes.
This validates table values and coverage, not biological provenance or accuracy.

Verification: test_reference_preflight.py (2 tests),
test_shared_cohort_runner.py (1 test),
test_target_recovery_reference_comparison.py (9 tests), and
test_metaphlan_reference_runner.py all pass. Shell syntax and git diff --check
pass. The new tests exercise the read-only reference helper and statically
check runner ordering; they are not a full-run end-to-end model test.

Items 1 and 2 above have local fixes, pending real-cluster verification.
Items 3–6 remain. Existing passed cluster preflights predate these stronger
reference gates and do not prove genome-size input readiness. No commit, push
or cluster submission was performed in this implementation step.
