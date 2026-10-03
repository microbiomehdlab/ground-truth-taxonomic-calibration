# Explaining missing adenoma significance with controlled spikes

Status: user-agreed scientific priority, 3 October 2026. This is a transparent
post-pilot clarification. It supersedes conflicting priorities in older plans.
It changes the purpose and reporting order; it does not rewrite completed
analyses, their settings, or their evidence status.

## Question and evidence chain

To what extent can the taxonomic profiling step help explain why CRC-associated
taxa do not consistently reach statistical significance in adenomas?

Use known read additions to measure what signal enters profiling, what signal
emerges, and whether it reaches significance under a fixed established DA
workflow. Relate those experimental requirements to the profiler-observed
adenoma distributions. CRC backgrounds remain mandatory secondary evidence.
The intended endpoint is an explanation of limitations, not a new biomarker
list, a ranking of profilers, or development of a statistical method.

An adenoma association can be absent, smaller, heterogeneous or poorly measured.
Spiking does not establish that a genuine adenoma disease effect exists.
Successful recovery is informative and must remain in the explanation.
Low absolute abundance, small group difference, low prevalence and small sample
size are distinct explanations; do not merge them into “low abundance.”

## What each component contributes

| Component | Question | Existing evidence | Remaining deliverable |
|---|---|---|---|
| Clinical reference / DA1 | What is the observed adenoma pattern relative to CRC and controls? | Sealed original profiles, 511-person descriptive reference; twelve clinical pilot contexts | Compact ten-target clinical context table with effect, uncertainty, p/q, prevalence, group sizes and observed distributions |
| Recovery | Does the profiler recover the known addition, and how consistently? | Three-cohort Bracken all-input and MetaPhlAn genome-equivalent endpoints and reports | Expose sample-specific recovery and native changes for the exact selected experiment participants |
| Paired response / DA2 | Is the added signal distinguishable within the same person? | Paired pilot and diagnostics; original/spiked profiles | Dose-response summaries for every target; keep paired significance supplementary to magnitude/consistency |
| Between-person response / DA3 | Does a known addition remain significant amid real between-person variation? | Deterministic 100-allocation plan, existing dose/n/heterogeneity arms, batched runner saving parametric and permutation results | Verify the completed collection, then report standard-workflow significance alongside recovery and background variation |
| Clinical interpretation | Are experimental limitations relevant to adenoma abundance ranges? | Native clinical distribution tables and dose/reference audits | Profiler-specific displays relating native baseline/spiked ranges to clinical distributions, with explicit scale limits |

No upstream rerun is required. The user's last messages do not establish the
live completion state of the large DA3 batch; inspect its report when supplied.

## Statistical role, fixed before the next synthesis

Use the established MaAsLin2 workflow as the clinical DA reference: the existing
pinned backend, native fractions, frozen transformation and full-family BH.
Keep the existing age/sex adjustment for actual disease contrasts and report
covariate eligibility. Do not choose settings to obtain adenoma nonsignificance.

The current DA3 worker already writes `parametric_p` and `parametric_bh_q`
alongside permutation results. These come from the validated fast group-only
equal-variance LM calculation, not a fresh call to MaAsLin2 for each context.
Describe that implementation accurately. Model/backend equivalence checks are
already available. Use these existing outputs for the standard-workflow view.
Do not relabel maximum-statistic FWER p-values as BH q-values.

Known null inflation under some small-n/unequal-variance settings remains
reportable evidence about the DA procedure. Standard usage does not make those
flags disappear. Retain existing permutation/HC3 comparisons to show whether an
interpretation depends on inference, without opening another method-development
campaign. Stop further wild-bootstrap expansion for this objective.

Native abundance is the DA input. Bracken all-input recovery and MetaPhlAn
genome-equivalent expected recovery are separate measurement quantities.
The first-20-feature wild-bootstrap output is an engineering subset; it supplies
neither final target inference nor a validated full-family clinical result.
Paired mixed-model failures remain documented; they do not require an unpaired
substitute for same-person measurements.

## Immediate implementation: link recovery to the exact DA3 experiment

Implementation update, 3 October: `link_da3_recovery.py` and its snapshot-based
launcher now implement this join and pass end-to-end local fixtures. See
[DA3_RECOVERY_LINK.md](DA3_RECOVERY_LINK.md) for exact inputs, output fields and
the cluster command. Real-input execution remains pending. Parametric SE is
reconstructed from selected native values; confidence intervals are not added.

The current `audit_recovery_vs_discovery.py` explicitly uses
`BACKGROUND_POOL_NOT_SELECTED_CASES` for DA3 recovery. Its 466 pilot rows are a
useful descriptive precursor, but cannot assign recovery to the selected people.
The large DA3 collector retains target statistics and context summaries; the
immutable plan retains the person/dose assignments needed to make the join.

Build two related source tables, with a shared interface across all cohorts:

1. **Person/target response table.** Key by cohort, profiler, biological sample,
   population, implanted target or community label, physical dose/profile identity
   and measured target. Include achieved total and target read fractions,
   original and spiked native abundance, native difference, profiler-appropriate
   expected added signal, recovered added signal, residual and recovery ratio.
   Carry reference type, input checksums and exclusions. Preserve negative/zero
   responses. A zero added dose has no defined recovery ratio.
2. **Context/target explanation table.** Key by context ID and target label.
   Join the plan's exact selected people and per-person doses to that response
   table and the collected DA statistics. Include cohort, background, profiler,
   n per group, arm, allocation, planned and achieved exposure, family size,
   transformed group effect, raw p, BH q and separately named permutation values.
   Summarize baseline abundance/prevalence/variation in both selected groups,
   native paired changes among exposed cases, recovery median/IQR and count of
   missing/zero/negative responses. Retain separate counts for exposed and all
   cases. Unexposed cases are not “failed recoveries.”

Join rules: exact frozen target aliases; sample ID includes cohort; independent
and community profiles cannot be interchanged; nominal total dose is not achieved
per-target dose. Variable-dose contexts join each person's assigned dose.
Missing/corrupt source evidence is an error, not a zero recovery or nonsignificant
result. Validate uniqueness and all intended participants/targets before summaries.
Never replace a missing individual endpoint by a pool median.

The compact DA3 target table does not retain parametric SE/CI or the explicit
parametric non-estimability status. Do not manufacture them from permutation
fields. Where required for explanation, reconstruct the fixed standard model
from the saved allocation and native profiles, with equivalence checks against
the retained p/q and unchanged full family. This needs no new permutations.

## Interpretation and attribution

Examine continuous recovery, variability and effect sizes before assigning
descriptive labels. Do not invent good/bad recovery cutoffs after inspecting
significance. Summarize within matched cohort/background/target/dose/n/arm
contexts; pooled correlations across doses are not a causal explanation.

- Weak/inconsistent measured additions accompanying nonsignificance are evidence
  consistent with a profiling contribution.
- Consistent additions with nonsignificance point to remaining effect-size,
  between-person variation, sample-size or multiplicity limitations.
- A small raw p but nonsignificant full-family BH q documents the effect of
  multiplicity for that contrast. A raw p is not a measure of recovery accuracy.
- Consistent recovery and significance limit how far profiling can explain the
  clinical observation under the tested conditions.

These comparisons do not isolate a causal fraction “due to profiling.” A stronger
claim would need an explicitly justified expected-response comparator holding
people, model, family and scale fixed. The existing counterfactual code is not
automatically validated for the current inputs. In particular, all-input Bracken
recovery cannot simply replace native classified-composition DA values.
For now report measured response and statistical detectability together, with
attribution bounded by these limitations.

## Execution order and completion criteria

1. **Reuse and verify DA3 collection.** Read its status and checksum manifests;
   require coverage of the planned 120,000 contexts before reporting full-grid
   frequencies. Preserve missing tasks separately. No fresh statistical pilot.
2. **Implement the two-table join above.** Read sealed endpoint tables and the
   immutable DA3 plan; verify person identity, scale and dose with small meaningful
   fixtures, including variable-dose and partial-exposure contexts.
3. **Make the explanation summaries.** Dose/n significance frequencies under
   standard LM/BH, recovery distributions and selected-person background
   variation, plus existing statistical sensitivities. All ten taxa stay visible.
   Repeated allocations measure conditional frequency, not additional people.
4. **Connect to clinical context.** Use the existing descriptive reference and
   fixed DA1 outputs. Display measured native ranges within each profiler.
   Never equate inserted read fraction with MetaPhlAn native cellular-like
   fraction, nor spike dose with a known clinical disease effect.
5. **Render a coherent figure set and write results.** Clinical context;
   added-versus-recovered signal; recovery-to-significance dose/n panels;
   clinically interpreted synthesis. CRC and sensitivity detail can be
   supplementary, but are retained regardless of direction.

Acceptance is complete, interpretable evidence, including good recovery and
successful detection. Neither proving profiler failure nor reproducing a desired
adenoma p-value is an acceptance criterion.

## Evidence already checked

The local wild-bootstrap REPORT bundle
`clinical_wild_checks_20261002T235941Z_REPORT` verified its checksums. It contains
108 summary rows, 200 simulations per context/scenario/method and a 20-feature
engineering family. Unequal-variance any-discovery ranges: ordinary LM 4.5–17.5%,
HC3 3–7%, wild bootstrap 4–7.5%. Eight ordinary-LM contexts had unadjusted
95% Wilson lower limits above 5%; none did for HC3 or wild bootstrap. Paired
profiler contexts can share the same simulated design and are not independent
validation replications. These diagnostics remain supplementary and do not
establish full-family clinical validity or biological absence.
