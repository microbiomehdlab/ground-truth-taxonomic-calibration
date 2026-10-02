# Proposed production matrix — 2026-10-02

Status: concrete proposal for user approval, NOT a production launch approval.
This is a post-pilot amendment, not preregistration. Once approved, reconcile
the older statistical protocol and implementation plan with this document;
do not silently treat their conflicting MaAsLin2-primary instructions as current.

## Purpose and evidence

Separate native measurement, quantitative recovery, paired technical
detectability, actual clinical association and artificial between-person
discoverability. Adenoma is the main background; CRC is mandatory secondary.
No design is selected to obtain a preferred direction of significance.

Clinical reference: `clinical_target_reference_20261002T172936Z`, 511 people,
10,220 native target rows, 180 target/cohort/profiler/condition summaries.
Downloaded checksums verified. Reported prevalence is not biological prevalence.
For example Bracken Fnuc CRC positivity was approximately 35–38%, versus
2–18% in adenomas; MetaPhlAn positivity was lower. Do not use those estimates
as the true fraction of people carrying a taxon or as fixed simulation parameters.

Robustness evidence and limitations: `ROBUSTNESS_REVIEW_20261002.md`.
Permutation maximum-statistic inference was reassuring under tested global
exchangeable nulls. Arbitrary partial-null strong error control and covariate
adjusted clinical inference have not been established by those experiments.

## Shared contract

- Cohorts: Yachida 201, Feng 154, Zeller 156; both profilers.
- Frozen ten target identities and exact aliases; no post-result substitutions.
- DA always uses native species relative abundances, not corrected recovery.
- Keep genome-equivalent MetaPhlAn recovery and all-input Bracken recovery
  separate from native DA; retain historical raw/reference comparisons.
- Transform: log2(1 + native_fraction / 1e-8); pseudocount sensitivities
  1e-9 and 1e-7 as already evaluated, clearly labelled.
- One immutable plan, same person allocations and exposure assignments across
  profilers, all doses and methods. No subject occupies both artificial groups.
- Full-family definitions frozen from the appropriate baseline pool using the
  established prevalence policy plus all ten targets. No spike-dependent filtering.
- Explicit status for non-estimable, absent, tied and zero-change endpoints.
  Do not silently remove tests, failed tasks or targets from denominators.

## DA1 — actual clinical associations

12 main contexts: 3 cohorts × 2 profilers × 2 contrasts
(Adenoma–Control and CRC–Control), using all eligible original people.
Primary candidate: pinned MaAsLin2 fixed clinical group + age + sex, no artificial
pairing. Keep current eligibility/missing-covariate ledger; no imputation or
extra adjustment variables added automatically. BMI/batch sensitivity requires
an explicit completeness/confounding review and decision.

Report both effect directions, uncertainty, raw p, full-family BH q, baseline
prevalence, sample counts and exclusions. Existing group-only DA3 validation
does not authorize DA1; validate the clinical design matrix and actual backend.
Compare effect/sign agreement across cohorts and tools, not only significant lists.
No interpretation that an adenoma nonsignificant result proves biological absence.

## DA2 — technical within-person experiments

Community: all people, separately by clinical condition; all seven existing
nonzero doses. Independent: all ten targets, six existing nonzero doses,
ten people per condition per cohort. Each comparison matches original and
spiked profiles by biological ID. Independent spikes cannot support sample
sizes above ten per condition.

Main outputs: distributions of native changes; positive/negative/tied fractions;
tool-appropriate quantitative recovery and reference uncertainty separately.
Direction sign test is a candidate inferential output, explicitly testing
equal sign probability among nonzero changes. Mean-change t and paired
MaAsLin2 results are comparisons with documented instability/skewness caveats.
Do not describe sign testing as evidence about mean change or recovery accuracy.
Retain full-family BH outputs; ten-pair discreteness/power limitations are explicit.

## DA3 — artificial independent-person discovery

All arms run in both adenoma and CRC backgrounds in every cohort and tool.
Use the existing community profiles and native baseline profiles only.
Artificial groups are assembled within a background: they are NOT observed
adenoma versus healthy-control comparisons. DA1 addresses that clinical contrast.

| Arm | Fraction of artificial cases receiving spike | Dose rule | Role |
|---|---:|---|---|
| U: uniform | 100% | One common existing dose per context | Favourable controlled benchmark |
| P: partial exposure | 25%, 50%, 75% | One common existing dose among exposed cases | Heterogeneity sensitivity |
| V: variable dose | 100% | Seeded equal-probability choice from existing low-dose grid per person | Heterogeneity sensitivity |
| PV: partial + variable | 50% | Same variable-dose rule among exposed cases | Combined heterogeneity sensitivity |
| N: no added signal | 0% | Native baseline in both groups | Null diagnostic |

An unexposed artificial case uses that person's original profile. An exposed
case uses that person's existing community-spiked profile. Every artificial
control uses its own original profile. Never interpolate or fabricate an
observed profile. Assignment must not depend on baseline abundance, measured
recovery or discovery results. All ten community targets are co-exposed;
partial exposure is community-level, not independent per-target exposure.

### Dose axis

Existing nominal community TOTAL read fractions:
0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1.
Per-target nominal amount is total / 10 for the equal-weight panel.
Thus total 0.0001 means 0.01% total and 0.001% per target.

U uses all seven doses. P uses the five lower doses (0.0001 through 0.01).
V/PV propose the fixed grid {0.0001, 0.0005, 0.001, 0.005, 0.01}.
The upper doses 0.05 and 0.1 are stress/positive-reference conditions, not
assertions of realistic clinical abundances. Variable-dose grid is a broad
sensitivity design, NOT a fitted disease model. It spans amounts rather than
being chosen to match significance. Publish actual achieved per-target amounts,
read additions and exposure counts per person/context.

Native clinical fractions are not interchangeable with inserted read fractions,
especially for MetaPhlAn. Do not use ratios of nominal read dose to native
clinical difference as accuracy or biological fold-change estimates.

### People, repeats and seeds

Community n/group = 5, 10, 15, 20; all backgrounds have at least 40 people.
Five/group is a deliberately small-n resolution diagnostic, not the only
main inferential result. Independent Fnuc/all-target ten-person subsets remain
DA2; a separate subset DA3 analysis, if retained, is exploratory n=5/group only.

Propose retaining 1,000 seeded allocation repetitions, as in the existing plan.
Draw disjoint groups and use nested n where existing allocation policy permits.
Exposure subset is a seeded ranking of case IDs, shared across dose/tool.
Use nearest-integer p*n with half-up rounding, recording achieved count and
fraction (e.g. 25% cannot be exact for n=5). Nested size need not imply nested
exposure sets; disclose that distinction. Variable doses are deterministic
functions of seed/cohort/background/allocation/sample ID, not worker order.

Discovery frequency is conditional on the finite observed cohort and artificial
allocation scheme. Repeated reuse does not create new biological people.
Monte Carlo intervals are simulation precision, not population confidence intervals.

### Inference and multiplicity — explicit proposal

Leading main candidate: exact enumeration at 5/group; Monte Carlo permutation
at larger n, same transformed contrast and joint person-label permutation
across species. Full eligible-family maximum-statistic adjusted p <= 0.05;
positive direction required for an enrichment discovery. This is FWER, not BH q.
Preserve full-family parametric MaAsLin2-equivalent group-only BH as a comparison.

Keep 9,999 Monte Carlo permutations as the validated starting budget;
report precision near the decision boundary. Do not interpret low Monte Carlo
BH discovery as adequate power: its p resolution can preclude isolated BH hits.
Predefined ten-target-panel permutation BH may be an explicitly separate
sensitivity, never substituted silently for full-family discovery. Any proposed
primary target-panel switch needs approval and dedicated calibration validation.

Global-null validation is NOT proof of strong FWER under partial signals.
Before using full-family discoveries as definitive, include partial-null synthetic
checks or justify the needed subset-pivotality conditions. Report non-target
changes as off-target responses; compositional changes are not automatically
false positives. Do not claim general clinical FDR control from these diagnostics.

## Workload and parallel execution

U: 7 contexts per allocation/n/stratum; P: 3 × 5 = 15;
V: 1; PV: 1; N: 1. Total 25 × 4 n × 12 strata × 1,000 allocations
= 1,200,000 conceptual DA3 contexts before sensitivity methods.
This is a workload estimate, NOT a reasonable one-job submission.

Build the immutable manifest first and measure representative timing before
authorizing that volume. Vectorize across species, shard contexts independently,
reuse identical validated null results where inputs/policy hashes match, and
avoid refitting identical group-only models for each target. Start with a bundled
canary covering every arm, background, profiler and sample size; then measure
runtime/storage and agree resource budget. A reduced repeat count or staged
heterogeneity expansion is a transparent design amendment, not a hidden shortcut.

Use bounded Slurm arrays (initial proposed concurrency 30 subject to allocation),
single-thread numeric libraries unless explicitly budgeted, source/image/input
hashes, independent output directories, deterministic retries and after-any
collector. Failed/unfinished shards remain failures, not biological non-discoveries.
Summary-only production outputs plus selected diagnostics avoid saving millions
of redundant full model objects. Retain raw p/adjusted values and exposure ledger
needed for reproducibility and selected figure sources.

## Freeze checklist / approval gates

1. Approve U/P/V/PV roles, grids, rounding and repeat/runtime budget.
2. Approve DA2 magnitude versus sign estimands and DA3 FWER/BH separation.
3. Validate DA1 covariate models and partial-null DA3 inference.
4. Implement immutable production plan + exposure assembly + reusable worker
   + fail-aware collector; test all arms together and benchmark costs.
5. Update current protocol/runbook without deleting historical decisions.
6. Only then launch production. This document does not assert those steps
   are already implemented or authorize new upstream generation.
