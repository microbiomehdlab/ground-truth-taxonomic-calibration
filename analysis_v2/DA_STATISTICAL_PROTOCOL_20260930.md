# Downstream statistical protocol: implementation baseline

**Amended 2026-10-01:** Read `DA_PROTOCOL_AMENDMENT_20261001.md` first.
MaAsLin2 is the agreed primary tool for all three DA analyses; custom tests
below are sensitivity methods. Historical text is retained, not current primary
model policy. Transformation/family implementation choices still require review.

Date: 30 September 2026. Scope: local protocol design, no analysis execution.

This resolves the proposed DA3 defaults in `DA_IMPLEMENTATION_PLAN_20260930.md`.
Use this document for statistical choices and that document for engineering,
outputs and tests. Defaults below are selected for implementation and validation,
NOT a claim that code, input feasibility or empirical calibration have passed.
Final execution freeze requires the gates below. Amendments must be dated and
justified; this redesign is informed by previous development, not preregistered.

## 1. Decision summary

| Item | Implementation decision |
|---|---|
| Main DA3 intervention | Community, every sealed dose |
| Backgrounds | Adenoma primary; CRC mandatory secondary |
| Independent comparison | Individual spikes and matched community subset, five versus five |
| Community group sizes | 5, 10, 15, 20 per group wherever feasible; no extrapolation |
| DA3 model | OLS log2(fraction + 1e-8) ~ artificial_case, HC3 covariance |
| Inference | Two-sided t reference, residual df=2n-2; 95% pointwise intervals |
| Discovery | Positive coefficient AND BH q<=0.05; q<=0.10 sensitivity |
| Feature universe | Background-specific full eligible baseline pool, prevalence >=10%, plus all ten panel targets |
| Multiplicity | Full frozen feature family separately per model context |
| Allocation repetitions | 1,000 independent seeded draws for full-pool community |
| Matched ten-sample subsets | Enumerate all 252 labelled five-case/five-control assignments |
| Uncertainty | Conditional allocation summaries; Monte Carlo error only, no population CIs |
| Primary surfaces | All tested doses/n; no required 80%/90% cutoff |
| Ideal-reference arm | Deferred, not a prerequisite for the main DA3 result |
| Full-run status | HOLD until implementation/input/null-validation gates pass |

HC3 and a t reference are an approximate small-sample procedure, not a guarantee
of calibration for sparse taxa. The null-validation gate is substantive. If it
fails, amend the protocol before definitive spiked analyses; do not method-shop
using target recovery. No new MaAsLin requirement is introduced by this protocol.

## 2. Roles of the three DA analyses

DA1: native baseline disease association, existing age/sex model and HC3 policy,
CRC/Control and Adenoma/Control contrasts. Existing BMI sensitivity remains.
Neither native significance nor absence establishes causal biological truth.

DA2: paired original/spiked samples, condition-stratified primary summaries;
existing two-sided one-sample test of paired log2 changes. All ten independent
samples in a condition can contribute paired observations. Do not apply the
five-per-group DA3 limit to DA2. See newly discovered code-policy gaps below.

DA3: distinct artificial groups drawn WITHIN adenoma or WITHIN CRC. Case/control
means assigned intervention state, not clinical disease/control classification.
Use cases' spiked profiles and different controls' baseline profiles. Full-pool
community and matched-subset experiments are separate, clearly labelled outputs.

## 3. DA3 estimand, model and numerical rules

For feature j, y_ij=log2(a_ij+1e-8), with native profiler abundance a as a fraction.
Fit y=intercept+beta*case. Beta is a difference in mean transformed abundances,
not arithmetic mean abundance, absolute organism counts or a clinical odds ratio.
Use no primary covariate adjustment because assignment is randomized. Retain
clinical metadata and report balance; never rerandomize to obtain significance.

Use sandwich HC3 and t with residual df=2n-2, matching the inference convention
already used by the native disease implementation. Cross-check against lm plus
the reference sandwich implementation. Store df, standard error and pointwise
intervals. These intervals are not simultaneous post-selection confidence bands.

Validate finite nonnegative fractions in the permitted native scale. Do not
silently species-close, multiply percent data twice, impute missing profiles or
clip unexpected results. Species absent from a COMPLETE profile are zeros only
under the validated parser contract. An unreadable/missing file is fatal.

Zero/near-zero residual variance needs deterministic handling: use a documented
scale-aware tolerance and test it on fixtures (candidate tolerance 100*machine
epsilon*max(1,max(abs(y))) for residual magnitude). A fit with residuals entirely
below that threshold is non-estimable for this protocol: retain beta, report SE/
interval/raw p as NA, p_for_BH=1, q bookkeeping value, status and reason. Never
manufacture an extremely small p-value solely because residual variance is zero.
All-zero features remain intended tests and never count as discoveries.

Unexpected rank failure, negative/nonfinite covariance, invalid values or an
unexplained numerical exception stop the context. Do not convert software bugs
into negative results. Ordinary estimable tests use two-sided p; do not halve
p-values after selecting the observed effect direction.

## 4. Feature family and BH

Freeze once per cohort/profiler/background from the full eligible COMMUNITY
baseline pool: species with nonzero baseline abundance in >=10% of its samples,
union all ten panel targets. Store prevalence numerators/denominators and the
exact frozen taxon mapping. No post-spike abundance is used to select features.

Use this SAME family for full-pool community, individual and matched-community
comparisons at every dose/n/allocation/null arm. This avoids changing the BH
family just because the matched subset has ten people. It is conditional on
these observed cohorts, not prospective training-only feature selection.

In individual experiments all ten panel species remain in the family, but only
the individually implanted species is an intended positive. In community all
ten are intended positives. Keep all other eligible native features in BH.

BH family key: cohort/background/profiler/pool/intervention/implanted-target-if-
individual/dose/n/allocation/arm. Fit one community contrast and extract ten
targets, rather than refit for each target. Degenerate tests contribute
p_for_BH=1 so they do not shrink multiplicity. Raw p remains NA and estimability
is explicit. Do not present bookkeeping q as a successful inferential estimate.

q<=0.05 and positive beta defines target discovery; q<=0.10 is sensitivity.
Negative significant results have a separate flag. BH is per contrast, NOT
global FDR over all doses, taxa surfaces, backgrounds and repetitions. Claims
about searching many cells require a separate inferential analysis; the main
surfaces describe prespecified algorithm performance without such claims.

## 5. Allocation and dose design

For full pools independently draw a random permutation per cohort/background/
repetition. If m is maximum feasible planned n, first m entries form ordered
case pool, next m control pool. First n in each forms each nested experiment.
Use the same assignment across profilers and every dose. No person can occur
in both groups; source profile versions do not create additional people.

Use n=5/10/15/20 only where 2n<=eligible N. Cross-cohort matched-n summaries
use common support; retain cohort-specific valid cells separately. Larger
cohort-specific n is deferred, not needed for the primary grid. Confirm CRC
counts and all dose availability from sealed inputs before manifest freeze.

Dose remains achieved added read fraction. Community total and per-target dose
must be separate fields. Do not assume equal weights mean exact total/10.
Use all available doses within cohort, common validated dose IDs for joint
comparisons, and report achieved ranges. Individual/community matching requires
overlapping per-target addition; exact matching is preferred. Any nonzero
tolerance requires a dated input-based justification before outcomes. If no
match exists report separate curves, not a fabricated matched comparison.

## 6. Finite individual experiment: enumerate, do not pretend 1,000 unique splits

Ten samples yield choose(10,5)=252 labelled assignments. Cases and controls
have different profile states, so complementary assignments are distinct.
Enumerate all 252 deterministically for individual and matched community,
with identical assignments across targets/profilers/doses. This is evaluation
over possible study allocations, NOT a permutation p-value test.

Report exact conditional frequencies over those assignments, with no Monte
Carlo intervals. Enumeration does not create population-level uncertainty or
increase biological N. Community full-pool n=5 remains a different experiment.
If subset N is not ten, stop and revisit the design; do not silently substitute.

## 7. Reproducibility estimands and denominator rules

For each target/context retain D_r in {0,1}: valid positive discovery, or not.
All planned valid-input allocations enter the denominator, including declared
non-estimable target tests as no discoveries. Input/software failures are fatal,
not excluded draws. Report estimability rates alongside discovery rates.

Within cohort, profiler agreement uses identical assignments: mean(D_K*D_M)
for both; the other three patterns derived from the same paired indicators.
Do NOT multiply profiler marginal rates because assignments are shared.

Across cohorts, allocations are independent by design. For each target primary
joint frequency is the PRODUCT of cohort marginal conditional discovery rates,
equivalent to averaging all cross-cohort allocation combinations. This avoids
an arbitrary pairing of replicate numbers. For example P(both)=f_A*f_B and
P(A only)=f_A*(1-f_B). All eight patterns use the corresponding products.
This product is justified by independent allocation, not independent biology
or a claim that real-world cohort outcomes have no shared biases.

All-three joint rate is f_A*f_B*f_C. Average target-specific joint rates over
the fixed ten targets for the panel summary; do not multiply panel mean rates.
For a single individual intervention the denominator is its intended target.

Directional panel replication A->B is sum_t(f_A,t*f_B,t)/sum_t(f_A,t);
zero denominator is NA. This is discovery-weighted conditional replication,
not the mean of replicate-level ratios or an additional inferential validation
test. Report numerator, denominator and both directions. For same-cohort
profiler replication use empirical joint rates in the numerator instead.

Native DA1 cross-cohort sign/effect agreement is separate; synthetic replication
does not declare native disease biomarkers true or false. CRC-minus-adenoma
frequency differences are descriptive, not paired subject contrasts or causal
effects of disease. Ten chosen taxa cannot support universal species claims.

## 8. Conditional uncertainty and stopping

Community R=1,000 draws with replacement from the space of valid allocations,
independent per repetition. Reused people induce no additional biological N;
independent random draws allow CONDITIONAL Monte Carlo precision estimates.
For a single-target rate report MCSE=sqrt(f*(1-f)/R), maximum about 0.016 at
R=1,000, and optional Wilson Monte Carlo interval clearly labelled. Shared
draws correlate estimates across profilers/doses/n; do not use unpaired SEs
for their differences. Compute differences at the allocation level for MCSE.

For cross-cohort products use cohort-wise resampling of allocation indices to
quantify Monte Carlo precision if intervals are wanted; retain within-cohort
pairing across profiler/dose/n. These are NOT biological bootstraps. Enumeration
has no Monte Carlo error. Population CIs and outer biological bootstrap are
out of scope for the first implementation; do not invent them in figure code.

No efficacy-based early stopping. If compute is insufficient, amend R globally
before definitive outcomes. Do not increase R selectively for attractive cells.
Full frequency surfaces are primary; threshold contours are deferred.

## 9. Null validation and pilot gates

For every planned allocation/n/background/profiler/pool evaluate baseline in
BOTH groups, same feature family/model. Null is dose-independent and fitted
once, not duplicated across dose denominators. Report any q<=0.05 discovery
in EITHER direction, number of discoveries, target positive-call frequency,
non-estimability and failures. Under this global randomized-label null,
any-discovery frequency equals empirical FDR; per-feature rejection fraction
is not FDR. This does not establish FDR for compositional spiked alternatives.

Stage A: synthetic fixtures, zero/constant/positive/noisy cases, reference model
agreement, BH hand checks, identities, allocation and dose checks.
Stage B: 20-allocation real-input runtime/schema smoke pilot (not enough to
judge calibration); no biological conclusions.
Stage C: full R=1,000 baseline-null evaluation and all252 subset null fits,
before definitive spiked fits. Review all contexts. Flag any null any-discovery
rate >0.075 at q=.05 or failures; also flag a 95% Monte Carlo lower bound above
.05 where applicable. These are engineering review flags, not formal global
hypothesis tests. Enumeration uses exact conditional rate, no MC bound.

Flags block definitive launch pending documented statistical review; lack of
a flag is not proof of calibration or power. Record conservative behavior too.
Any method revision follows null diagnostics/theory and reruns gates, never
selects the method giving the best spike story. Correct software failures first.

In spiked arms significant non-targets are NOT automatically false positives:
relative dilution can be a genuine intervention effect. The baseline null is
the required negative control, not a license for off-target removal.

## 10. Existing-code findings: do not silently describe them as fixed

Inspected `scripts/fit_paired_biomarker_models.R` (canonical, not just its fast
fork). It calculates prevalence over baseline PLUS all dose matrices. Thus
the historical description of baseline-only paired filtering is incorrect.
It also assigns machine-minimum p to constant nonzero paired changes in
safe_test(). Neither behavior was changed during this review.

Required DA2 review: approve a dated amendment to baseline-only filtering and
conservative non-estimability handling, or explicitly retain and justify the
historical approach. Do not publish DA2 as baseline-filtered until implemented
and tested. DA3 baseline-only rules above are new, not proof of DA2 compliance.
Do not treat a constant nonzero response as an ordinary valid t-test silently.

Native disease code uses baseline-derived features and explicit HC3/t inference.
Its family/error handling must still be checked under real data. Existing
ANALYSIS_POLICY.tsv includes historical calibration/ideal-reference roles;
it was NOT rewritten here. A versioned policy migration must distinguish
retired scope from required main analyses without invalidating historical runs.

## 11. Implementation handoff and final freeze checklist

Implementation may begin from these defaults; no full cluster runs yet.

- Inventory final cohort/background N and doses, target-dose match availability.
- Resolve DA2 amendment and old/new policy version mapping explicitly.
- Repair runner/reference/consumer blockers from the scope audit.
- Implement and test DA3 model, constant handling, fixed feature family, 252
  enumeration, nested random allocations and replication formulas above.
- Produce human-readable counts/costs/config and execute the three pilot gates.
- Record code/container/data hashes, seed scheme, final manifest and all gate
  reports; mark final-execution status only after review.
- User runs cluster commands. No cluster access, manuscript changes, commits,
  pushes or data uploads performed by this document update.

## Primary references checked for this decision record

- R stats t.test documentation (constant data and test behavior):
  https://www.stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html
- R stats p.adjust documentation (BH and padding family with p=1):
  https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
- sandwich package documentation (robust covariance implementations):
  https://sandwich.r-forge.r-project.org/articles/sandwich.html

These references document mechanisms; they do not validate the proposed model
for these microbiome data. Empirical and implementation gates remain required.
