# DA protocol amendment: MaAsLin2 primary

DA2 primary-tool decision below is superseded by
[the 2026-10-02 paired-inference amendment](DA2_PAIRED_INFERENCE_AMENDMENT_20261002.md).
DA1 and DA3 remain MaAsLin2 primary. The dated original decision is retained.

Status: agreed primary-tool change; remaining proposed settings require review.
No definitive DA launch, preregistration claim or implementation-completion claim.
This document supersedes the HC3-primary and one-sample-test-primary passages
of DA_STATISTICAL_PROTOCOL_20260930.md and DA_IMPLEMENTATION_PLAN_20260930.md.
Their allocation, dose, uncertainty and replication definitions remain unless
explicitly amended below. Historical code and ANALYSIS_POLICY.tsv are unchanged.

## 1. Agreed decisions

- MaAsLin2 is primary for DA1, DA2 and DA3; custom tests are sensitivity analyses.
- Use the pinned MaAsLin2 1.18.0 environment, with real image/version verification
  before fitting. No silent update or replacement with MaAsLin3.
- Analyse cohorts and profilers separately; compare results afterwards.
- Quantitative recovery reference adjustments do not modify DA input abundance.
- Community DA3 is main; adenoma background primary, CRC mandatory secondary.
- Individual spikes and matched community subsets provide complementary evidence.
- Keep full-family multiple testing; do not test only the implanted target panel.
- Shared inputs, allocation rules, outputs and commands across all three cohorts.

## 2. Primary models

| Analysis | Fit context | Fixed effects | Random effects |
| --- | --- | --- | --- |
| DA1 | Separate CRC/control and adenoma/control comparisons per cohort/profiler | clinical_group + age + sex | None |
| DA2 | Separate population/intervention/dose/clinical-background contrasts | spike_state | biological_sample_id intercept |
| DA3 | Separate allocation/background/intervention/dose/n/profiler contrasts | artificial_case | None |

DA1 control is the reference category; retain both effect directions, extract
the ten panel taxa alongside full native results. Adjustment availability must
be audited: no silent covariate dropping, metadata imputation or automatic
sample removal. A reviewed complete-case eligibility ledger is required before
family freeze. BMI sensitivity is separate and conditional on availability.
Batch adjustment is not added automatically; audit confounding and support
before deciding any cohort-specific sensitivity.

DA2 has one original and one spiked version per person per fit, not all doses
treated as independent samples. Use MaAsLin2 analysis_method=LM with sample ID
as random_effects to handle pairing. Ten individual samples per condition are
ten biological pairs. Community can use all eligible pairs in that condition.
Time-invariant age/sex are not required in the primary within-person fit.

DA3 has different people in each group and only one profile per selected person.
No artificial sample-ID random effect; randomized assignment supports primary
group-only fitting. Retain clinical metadata and balance, without rerandomizing.

## 3. Transformation: AGREED 2026-10-01; implementation validation pending

Recommend retaining y=log2(native_fraction + 1e-8), identical across contexts,
with native percentages converted exactly once. Pass a nonnegative shifted
version y-log2(1e-8) to MaAsLin2 with normalization=NONE and transform=NONE.
The fixed shift only changes the intercept, not group effects; validate this
with direct LM/mixed-model fixtures. This is an explicit non-default MaAsLin2
workflow, not its native LOG transform, and must be described accordingly.

Reason: the current upstream MaAsLin2 LOG implementation replaces zeros with
half the minimum positive abundance of each feature, then logs it. Re-estimating
that minimum in each resampled context can change zero handling across draws,
doses and n. Using built-in LOG instead is a valid alternative if explicitly
chosen and verified in pinned 1.18.0; do not mix the approaches unnoticed.
AST is another substantive alternative, not an automatic copy of another paper.

No TSS reclosure after feature selection: native outputs are already relative
abundances on profiler-specific denominators. No learned calibration, recovery
ratio input or post-spike-based pseudocount selection. Use standardize=FALSE
for deterministic units; binary/categorical effects have explicit references.
The 1e-8 pseudocount is a modelling choice, not a measured detection limit.

User agreed to this primary transformation with sensitivity checks and baseline
null validation. Do not claim it is uniquely correct or reviewer-proof.
Prespecified pseudocount sensitivities: 1e-9 and 1e-7, plus pinned MaAsLin2
built-in LOG as a separately labelled sensitivity. Keep the same eligible
samples/family/allocations and never choose the primary based on target outcomes.
Sensitivity scope and cost must be inventoried before submission.

## 4. Feature families, filtering and multiplicity: proposed implementation

Freeze native species families from baseline prevalence >=10%, union the ten
targets, independently for each cohort/profiler and analysis eligibility pool.
DA1: comparison-specific eligible disease/control baseline pool.
DA2: clinical-background full community baseline pool.
DA3: same full clinical-background baseline pool used in the existing plan.
Use the DA2/DA3 baseline family also for independent and matched community
subsets; record identities and numerator/denominator prevalence counts.
This is conditional baseline selection, not prospective training-only selection.

Disable context-dependent internal prevalence/abundance filtering after this
external freeze (min_prevalence=0, min_abundance=0). Record min_variance=0 and
test which constant features the pinned package omits. Do not assert these
settings prevent every internal omission. Preserve all family members in the
wrapper's final result ledger.

Primary BH is recalculated externally from the p-values for the intended GROUP
coefficient across the FULL frozen feature family of that context. MaAsLin2's
native q-values remain in a separate field because it can adjust across other
metadata coefficients and omit features. Never label wrapper q as native q.
Use p_for_BH=1 for expected non-estimable features, raw p=NA and explicit reason.
Unexpected fitting/numerical failures block the context, not shrink its family.

Primary target discovery requires valid positive group coefficient and wrapper
BH q<=0.05; q<=0.10 sensitivity. Native DA1 reports both directions; directional
target replication is separate. Two-sided testing; no post-hoc p halving.
All-zero targets remain intended positives in synthetic denominators but cannot
be discoveries. Significant non-targets under spikes are not automatically
false positives because relative dilution is a genuine intervention response.

## 5. Variance and fitting diagnostics

Do NOT describe MaAsLin2-primary inference as HC3. Primary p/SE come from its
LM or mixed-model backend. Determine exact df and CI extraction behavior in
pinned 1.18.0; output unsupported quantities as NA rather than invent them.
Save selected model objects for validation and explicit model/status records.

All-zero/constant response or scale-aware near-zero residual variance receives
non-estimable status and conservative family bookkeeping. In paired fits,
singular random-intercept fits and convergence warnings need classification:
hold context for review rather than silently switch to unpaired inference or
turn warning-driven omissions into zero discoveries. Freeze tolerances and
acceptable diagnostic policies using synthetic and baseline-null gates before
definitive spike fits, not by choosing attractive target outcomes.

Sensitivity: existing baseline disease HC3 fits, paired tests, and DA3 HC3 fits,
after their known feature/filter/constant handling defects are corrected and
tested. Match samples, families and transformation wherever scientifically
possible, and label any differences. Agreement is robustness, not truth.

## 6. Allocation and dose plan retained

Community n=5/10/15/20 per group where feasible, 1,000 seeded allocations,
nested n and identical assignments across profilers/doses/null arms.
Ten-sample independent subsets: enumerate all 252 labelled five/five allocations;
matched community uses the same people and allocations. No repeated person
across case/control groups within a fit. Dose matching uses target addition,
not total community dose. Keep unmatched curves separate, with honest support.

Observed grid: community total doses 0.01%, 0.05%, 0.1%, 0.5%, 1%, 5%, 10%;
independent doses 0.01%, 0.05%, 0.1%, 0.5%, 1%, 5%. Final inventory must verify
these across all samples. Community per-target fractions are reconstructed
integer contributions, not automatically exact total/10. No extrapolation.

Retain conditional allocation frequencies, not biological population CIs.
Within-cohort profiler joint discovery is empirical paired allocation frequency;
cross-cohort products use independent allocation draws and target-specific
marginals, not products of panel-averaged rates. Enumeration has no MC error.

## 7. Before cluster model submission

1. Transformation is agreed; user reviews proposed wrapper-family rules (section 4).
2. Inventory native species inputs, metadata, samples/doses and computational
   context counts; confirm all family/eligibility policies explicitly.
3. Implement shared MaAsLin2 wrapper with full audit/result mapping, pairing,
   package version checks and no silent sample intersection.
4. Test against direct LM/mixed models, zero/constant/singular examples, group
   directions, external BH, input identity, transform and percentage conversion.
5. Run authorized 20-allocation runtime/schema pilot; no scientific conclusions.
6. Baseline-null gates: original profiles in BOTH artificial groups, same
   allocations/model/family; full 1,000 draws and 252 subset allocations.
   Review any-discovery rates, conservative behavior, fit failures and MC error.
   Existing review flags (>0.075 any discovery at q=.05; MC lower bound >.05)
   remain engineering-review criteria, not formal global significance claims.
7. Freeze versioned policy/config only after review, then definitive spike fits.

## Sources checked 2026-10-01

https://github.com/biobakery/Maaslin2 (official manual: fixed/random effects,
normalization, transformations, results and q-values).
https://github.com/biobakery/Maaslin2/blob/master/R/utility_scripts.R
(current LOG zero-handling mechanism).
https://github.com/biobakery/Maaslin2/blob/master/R/Maaslin2.R
(current defaults; confirm installed pinned implementation before execution).
These sources describe software, not proof of small-n microbiome calibration.
