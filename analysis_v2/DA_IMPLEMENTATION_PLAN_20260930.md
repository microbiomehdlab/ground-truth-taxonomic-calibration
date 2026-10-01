# Differential-abundance implementation plan

**Amended 2026-10-01:** `DA_PROTOCOL_AMENDMENT_20261001.md` supersedes
HC3-primary model passages below. MaAsLin2 is primary; custom models are
sensitivity. Use the amended model/output contract for implementation.

Date: 30 September 2026. Status: agreed scientific scope; statistical details
marked PROPOSED await review and protocol freeze. No new analyses have run.

This document extends `BIOMARKER_REPRODUCIBILITY_PLAN.md` and supersedes its
description of CRC as merely contextual and individual experiments as optional.
Read it alongside `DOWNSTREAM_SCOPE_AUDIT_20260929.md` for existing code blockers.
It is a dated redesign informed by development results, not a preregistration.

## 1. Scientific structure

| Analysis | Comparison | Scope |
|---|---|---|
| DA1 native associations | CRC versus controls; adenomas versus controls | Original unspiked samples; all three cohorts and both profilers |
| DA2 paired artificial biomarkers | Spiked versus original version of the same sample | Individual and community spikes, available doses, clinical backgrounds separately |
| DA3 controlled discovery | Spiked artificial cases versus different unspiked artificial controls | Community main experiment; adenoma primary, CRC required secondary; individual matched comparison |

Original plan modules A1/A3/A4 correspond to DA1/DA2/DA3. Quantitative recovery
(A2) and explanatory summaries (A5) support them, not additional DA experiments.
The baseline randomized-label null is a required validation arm of DA3.

DA3 groups contain the SAME clinical background: adenoma versus adenoma or CRC
versus CRC. It does not test adenoma diagnosis. Clinical controls are not a
required third DA3 background. CRC is prespecified regardless of adenoma results.

DA2 uses paired observations and may use all ten individual samples per condition.
DA3 uses distinct samples and is limited to five per group in that subset.

## 2. Decisions and unresolved choices

AGREED:

- Three cohorts, both profilers; shared analysis interface and output schema.
- Community spikes are primary, with ALL available sealed community doses.
- Adenoma primary background and CRC secondary background.
- Vary sample size as well as dose; test every eligible feature, extract all targets.
- Same biological allocations across doses and profilers; no sample in both groups.
- Individual spikes provide a matched small-n comparison, target by target.
- No new upstream profiling. Filtering/calibration do not return as required stages.

PROPOSED, to freeze after feasibility/statistical review:

- Common n = 5, 10, 15, 20 per group, conditional on actual availability.
- DA3 log2(abundance fraction + 1e-8) ~ artificial case, HC3 inference.
- Baseline-frozen feature family, 10% prevalence plus intended targets.
- 1,000 allocations after a small runtime/null pilot.
- Specific degenerate-fit policy, inference tails, replication aggregation.
- Optional ideal-reference arm and outer biological uncertainty procedure.

Existing DA1/DA2 frozen policies are not silently changed by the DA3 proposal.
Before definitive execution record signed-off choices, rationale, code/config
hashes and amendments. Never choose settings by which yield favorable targets.

## 3. Evidence and prerequisites

User-reported common upstream seals cover 201 Yachida, 154 Feng, 156 Zeller
samples. Shared preflights passed: 3097753/3097754/3097755, with canonical rows
36,360/28,840/29,160. These establish input readiness, not successful DA fitting.

First resolve existing endpoint-directory conflicts, missing final genome-reference
coverage checks and stale quantitative consumers identified in the scope audit.
The existing shared runner is a starting point, not a completed DA3 implementation.
Keep upstream seals and the evidence package immutable. No agent cluster access:
provide commands for the user after local verification and explicit run approval.

## 4. Input inventory and eligibility

Read canonical manifests and sealed profile inventories, not a recursive scan
that might include quarantined or obsolete outputs. Required attributes:

- Cohort/study/sample identity, clinical condition, independent-subset membership.
- Profiler, baseline/community/individual state, individual implanted target.
- Nominal dose identifier, achieved total and per-target read fractions.
- Profile path/checksum, units, taxon mapping version and completeness.
- Availability of both profilers and intended doses, required reference inputs.

Validate uniqueness of profile keys, sample condition consistency and baseline
coverage. Missing/corrupt profiles are input errors, not zero abundance. Complete
zero measurements remain eligible. Any scientific exclusion needs a recorded
reason fixed before outcome analysis. Fail closed on unexplained missing inputs.

Confirm counts directly. Previously reported adenoma pools are Feng47/Zeller42/
Yachida67 before exclusions. Enumerate CRC counts; do not infer them from total N.
Verify ten independent samples per condition/cohort and their containment in the
community pool. The 30 independent samples are NOT 30 adenomas or 30 CRC samples.

Primary comparisons use a fixed eligible sample pool across doses/profilers.
If grids differ, retain all doses within cohort but compare cohorts only on
documented common support. Do not silently merge doses or drop incomplete samples.

## 5. Dose definitions — mandatory before fitting

Community dose is TOTAL added read fraction across ten taxa, not the amount of
each species. Derive achieved target fractions from recorded read contributions
and denominator definitions. Do not assume total/10: configured equal weights
alone do not establish identical achieved counts.

Create `dose_mapping.tsv` with cohort, experiment, target, nominal dose ID/value,
achieved total fraction, achieved target fraction, denominator definition and
source provenance. Preserve sample-specific values and report their ranges.
Use every sealed community dose; enumerate actual values from inputs rather
than hard-code a number of doses from presentation slides.

Use total community fraction on community axes, with per-target additions in
source tables and target-specific threshold reporting. Profiler-specific genome
equivalent expectations are separate from achieved read truth. Baseline null
is a separate arm, not another independently fitted spiked-dose observation.

Individual/community matched-dose comparisons require matching target addition,
not equal nominal total fractions. Freeze a justified tolerance before outcomes.
If no overlap exists, show separate curves and label the limitation; do not claim
an isolated community effect or interpolate beyond observed support. Record both
target and total dose even when target doses match.

## 6. Allocation algorithm

For each cohort/background/repetition generate a seeded random permutation of
eligible biological sample IDs. Create two disjoint ordered case/control lists.
For each n take the first n from each list: nested sample-size experiments.
Assignment must be symmetric and randomized, not accession/clinical order.

Cases use spiked profiles; controls use original profiles of DIFFERENT samples.
The other state of a selected sample is not an extra observation. No biological
sample can appear in both groups within an experiment. Reuse assignments across
profiler, dose and baseline-null arm. Cohorts/backgrounds use independent stable
seed namespaces. Repetition numbers across cohorts pair independent allocations,
not the same people. Avoid unstable language-runtime hash functions for seeds.

Proposed shared grid: n=5/10/15/20 per group. Reject cells with 2n>N; write the
feasibility matrix. Optional largest balanced n=floor(N/2) is secondary, not a
matched-n cohort comparison. Odd N leaves one out. All people may participate
across repetitions; not everyone must appear in every draw. Report participation.

Allocation ledger: allocation ID, cohort, background, repetition, n, sample ID,
group, within-group order, seed namespace, eligible-pool hash. Separate profile
selection ledger links allocation to experiment/dose/profiler/source profile.
Resume/cache keys include design and input hashes. Repeated use across draws is
permitted but does not increase biological sample size.

## 7. Community experiment

Within each cohort/background/n/dose/repetition/profiler fit ONE full-feature DA
analysis. Extract results for all ten intended targets from that fit. Do not fit
the same community contrast ten times or describe it as ten isolated interventions.
Targets already present at baseline remain intended positives.

Primary pool: all eligible community samples from that clinical background.
Report per-target discovery, effects, intervals, p/q, prevalence and family size.
Preserve unrecovered targets rather than showing only successful species.

## 8. Individual-spike and matched-community comparison

For each cohort/background use its ten independent samples. Assign five cases
and five controls, separately for each repeated allocation; reuse allocations
across targets/doses/profilers. Cases receive one target's individual spike.
Each implanted-target intervention needs its own full-feature DA fit.
Only that implanted target is a truth-positive in its individual experiment.

Run community comparisons restricted to exactly these same ten samples with
identical assignments and validated matched target doses. Reuse a community fit
to extract multiple targets when its context is identical. This is distinct
from the full-pool community n=5 analysis; label both clearly.

Differences may reflect compositional changes and measurement interactions,
not solely classifier interference. Unequal target doses cannot isolate the
effect of joint addition. Small-n results are complementary, not a substitute
for the full-pool dose/sample-size experiment.

## 9. Model interface and multiple testing

Context key: analysis ID, cohort, background, profiler, population, intervention,
individual target if applicable, dose, n, allocation, arm. DA1/DA2/DA3 share
output conventions but need not use the same statistical model.

Proposed DA3 response: log2(fraction+1e-8), coefficient artificial case minus
control, HC3 uncertainty. Validate percentage/fraction units. Randomization
supports no primary age/sex adjustment; record balance without redrawing for
favorable p-values. Review small-n null and fit behavior before model freeze.

Proposed feature universe: baseline 10% prevalence plus targets, fixed within
cohort/profiler/background from its eligible unspiked pool. Explicitly reconcile
this background-specific scope with existing whole-cohort policies before coding.
This is conditional/transductive evaluation, not prospective training-only filtering.
Keep family fixed across n/repetition/dose/arm; use the same family for individual
and matched-community comparisons. No filtering on spike response or significance.

BH is over the entire eligible family for each contrast, not only ten targets.
Preserve family size, identity and raw p. Optional target-only multiplicity is
a separately labelled diagnostic, never the primary replacement after results.

Retain coefficient, SE, interval, p, q, group counts, prevalence, family ID,
estimability and status for every feature. Freeze degenerate handling. Proposed
conservative bookkeeping: raw p/interval NA for declared non-estimable features,
p_for_BH=1 to retain family, explicit reason. Unexpected numerical failures stop
the stage. This rule needs sign-off; never silently reduce the family.

Primary target discovery: valid positive coefficient and q<=0.05; q<=0.10
sensitivity. Report negative significant calls separately. Freeze test tails;
do not select one-sided direction from outcomes. Complete all-zero targets stay
in intended denominators and cannot be discoveries. Distinguish missing inputs,
expected non-estimability and numerical errors. Conditional-on-estimability
summaries are secondary, not a way to hide target failures.

## 10. Baseline null and non-targets

Use the same allocation with original profiles in BOTH groups. Fit once per
allocation/n/profiler/background: no dose dependence. Reuse across dose displays
without counting duplicates in null denominators.

Report any-discovery frequency, discovery counts, target-specific false-call
frequencies and fit failures. Under the global randomized-label null, empirical
any-discovery frequency corresponds to FDR; fraction of features significant
does not. Results remain conditional on these cohorts and this procedure.

In spiked arms non-implanted taxa can genuinely decrease in relative abundance
through dilution. Do not label every significant non-target a false positive.
Formal non-target truth needs an explicit compositional expectation and null.
Do not introduce off-target filtering or learned output calibration here.

## 11. Replication endpoints

At matched n/dose/background record per-target profiler both/first-only/second-
only/neither; all cohort pairs; all eight three-cohort discovery patterns.
Always retain single-cohort frequencies and effect directions.

Joint recovery denominator: all intended targets in valid input contexts, with
failure flags visible. Directional replication denominator: discoveries in the
designated discovery cohort evaluable in validation. No discoveries means NA.
Specify pooled conditional replication versus mean replicate-level ratios;
they differ. Freeze primary aggregation and report actual denominators.

Use independently generated cohort draws; cohort pairs are not independent
biological replications of one another. Ten targets are a selected panel, not
a representative random sample of all species. Preserve profiler-specific
feature families rather than pretending q-values have identical multiplicity.

CRC runs use the same rules and common feasible grid. Present CRC separately,
then descriptive adenoma–CRC differences. These are different people, not paired
subjects; differences cannot establish disease as a causal driver of measurement
failure. Same added dose need not yield the same relative effect at baseline.

## 12. Repetition, uncertainty and computational plan

Begin with synthetic fixtures, then a proposed 10–20-allocation authorized
runtime/null pilot. Review resources, schemas, fit stability and null behavior,
not which taxa support the story. Freeze substantive changes with rationale.

Proposed definitive R=1,000. More repetitions improve Monte Carlo precision,
not biological N. Report conditional allocation frequencies. Do not present
binomial intervals over repetitions as population CIs. Outer biological
resampling is unresolved and needs its own design preventing an original
person from entering both groups. Identify duplicate small-panel allocations
without treating them as new biological data.

Let C=3, B=2, P=2, K=feasible n levels and D=community dose count. Community
cost is C*B*P*K*D*R full-family fits, NOT times ten targets. With K=4 this is
48*D*R. Baseline-null cost is 48*R. Individual cost sums D_target*R across
cohort/background/profiler/implanted target at n=5. Matched-subset community
and null add fits; deduplicate identical contexts. Each family fit evaluates
many features. Enumerate actual counts and benchmark before Slurm requests.

Reuse input matrices, not copied profiles per draw. Partition by cohort,
background, profiler and repetition range with atomic checkpoints. Merge only
complete nonoverlapping contexts; reject stale hashes. A failed chunk must not
alter seeds, planned repetitions or denominators. Resume must equal serial run.

## 13. Ideal-reference comparator and interpretation

Optional proposed arm: replace case profiles with expected post-spike profiles
on the appropriate profiler scale. Same allocation/model/family/BH. Validate
the complete compositional table and genome-reference coverage before use.
Do not call observed-minus-ideal an exact causal decomposition of nonreplication.
This comparator is not learned calibration.

Primary output is a frequency surface across TESTED dose/n, without forced
monotonicity. Any 80%/90% operational contour must be chosen beforehand and
labelled descriptive. No crossing means not reached/not established; do not
invent interpolation or extrapolate beyond sample/dose support. Native disease
log-fold effects cannot be directly converted into spike read fractions.
No manipulated-depth experiment exists, so infer no universal required depth.

## 14. Output specification

Names below are proposed interfaces, not claims that files already exist.

| Artifact | Content |
|---|---|
| design_config | Versioned agreed/frozen settings, seeds, model, tails, q, R and grid |
| input_inventory.tsv | Profile keys, cohort/population, provenance and hashes |
| eligibility.tsv | Counts, scientific exclusions, per-cell feasibility |
| dose_mapping.tsv | Nominal/achieved total and target fraction; denominator |
| feature_universe.tsv | Identities, baseline prevalence, target flag, family IDs |
| allocations.tsv | Biological group membership and nested n |
| selected_profiles.tsv | Exact input for every context |
| da_results partitions | Full feature coefficients, uncertainty, p/q, diagnostics |
| target_discoveries.tsv | Every intended target, positive/negative/failure states |
| replication_patterns.tsv | Profiler/cohort patterns, denominators |
| conditional_frequencies.tsv | Frequency, draw count, Monte Carlo labels |
| null_diagnostics.tsv | False calls under randomized baseline labels and failures |
| individual_community_comparison.tsv | Matched subset/allocation/target-dose flags |
| run manifest and checksums | Code/config/input/environment/output provenance |

One schema for all cohorts. Keep operational absolute paths separate from
shareable relative-path mappings where needed. Observe data-sharing terms;
this plan authorizes neither uploading nor publishing data automatically.

## 15. Figures and source data

1. Adenoma community target-by-target dose versus n discovery surfaces.
2. Matched-dose/n cohort/profiler joint recovery, with single-cohort frequencies.
3. CRC secondary with identical layout/scales and all targets visible.
4. Individual/community matched ten-sample comparison, n=5, dose-match flags.
5. Null diagnostics by n and background.
6. Supporting baseline/prevalence/error/variability descriptive summaries.

Use one profiler palette, taxon order, clinical labels, typography and unit
convention. Every panel has a rendering script, source TSV and caption stating
population, n per group, dose meaning, repetitions and uncertainty definition.
Never omit zero recovery or missing matched-dose support. Figure placement in
main versus supplement remains editorial, not a reason to skip planned results.

## 16. Mandatory tests

- Unique 2n samples; no cross-group overlap; nested n; stable seeds; no profile
  duplication inflating sample counts.
- Same allocation profiler/dose/arm; independent cohort namespaces; scheduling
  and resume cannot alter assignments.
- Missing profile fails; valid zero stays; duplicate keys fail.
- Dose total/target distinction; unequal dose cannot be marked matched; reject
  percent/fraction mismatch.
- Hand-calculated BH; fixed family; target-only analysis cannot pass as primary;
  negative effect cannot count as positive discovery.
- Exact degenerate policy; unexpected failures stop rather than shrink family.
- Synthetic positive/null model fixtures and correct paired/unpaired identity.
- Four/eight replication patterns; zero denominator NA; aggregation definitions.
- Matched independent/community sample sets; only individual implanted target
  counts as intended positive in that intervention.
- Null deduplication across dose; correct repetition denominator.
- Reject overlapping/missing chunks and stale hashes; resume equals serial.
- Three-cohort common schema/CLI and output-to-source provenance coverage.

## 17. Work packages and acceptance gates

| Order | Work | Completion criterion |
|---|---|---|
| 0 | Confirm actual doses, CRC counts and statistical choices | Reviewed feasibility matrix and dated frozen config |
| 1 | Fix known shared-runner/reference/consumer blockers | Real caller/helper and reference-scale tests pass |
| 2 | Inventory, dose/family ledgers, allocations | Fixtures plus human-readable design summary pass |
| 3 | DA3 full-family model/null and target extraction | Positive/null/degenerate/BH tests pass |
| 4 | Replication and matched individual comparison | Patterns/denominators/dose matching verified |
| 5 | Small authorized real-input pilot | Diagnostics and cost reviewed, protocol frozen |
| 6 | Definitive common cohort execution and aggregation | Complete contexts, diagnostics, stage provenance pass |
| 7 | Plots, methods and publication package | Source tables, consistent style, all planned results retained |

Stop before definitive execution for ambiguous dose meaning, missing inputs,
overlapping groups, unresolved feature policy, unexplained null/fit problems or
failing reference-scale tests. No upstream rerun is proposed. No implementation,
cluster submission, manuscript edit, commit or push has occurred for this update.
