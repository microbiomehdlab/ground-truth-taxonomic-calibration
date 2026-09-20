# Downstream v2 methods decision log

This tracked log records decisions that affect manuscript methods or
interpretation. Generated run directories preserve the corresponding inputs,
diagnostics, provenance, and checksums.

## 2026-09-21 — Yachida residual-genome-size audit completed and verified

**Result.** André resumed the sealed Yachida `DEVELOPMENT_ONLY` propagation at
`work/geff_propagation_dev_20260920T003139Z`. The runner reused every upstream
stage, quarantined the obsolete comparison intact, regenerated only the paired
reference comparison, and ran the residual-genome-size audit. It exited with
status 0, and every entry in the regenerated `run_checksums.sha256` verified.

The prespecified target-level regressions gave:

| Scope | Read-proportional sensitivity slope (95% target-bootstrap interval) | Genome-equivalent primary slope (95% target-bootstrap interval) |
|---|---:|---:|
| Independent | -0.9713 (-1.0669, -0.9181) | 0.0132 (-0.0754, 0.0612) |
| Community | -0.9995 (-1.1094, -0.9256) | 0.0068 (-0.0928, 0.0647) |
| Pooled | -0.9962 (-1.1019, -0.9247) | 0.0080 (-0.0903, 0.0642) |

All six fits used the ten implanted taxa as the statistical and bootstrap unit;
all 10,000 bootstrap replicates were valid. Read-reference fits had
`R² = 0.992–0.994`, whereas corrected fits had `R² = 0.008–0.037`. The
validation table records exactly one cohort (`yachida`), 15,870 MetaPhlAn and
15,870 Bracken observations, `bracken_identical = PASS`, measured spike-FASTA
genome sizes, no fitted genome-size constant, and no pseudocount. Paired
eligibility excluded 1,856 MetaPhlAn observations (11.7%) because at least one
arm had a nonpositive ratio; the per-arm reason counts overlap and must not be
summed.

**Interpretation.** Within Yachida, the near -1 genome-size relationship under
the read-proportional reference is almost completely removed by the
genome-equivalent reference. This strongly supports abundance-scale mismatch as
the dominant explanation for the earlier genome-size trend. It does not prove
that all residual taxon- or context-specific MetaPhlAn error is absent.

**Boundary.** This supersedes the earlier status statements below that the real
Yachida slope was unknown. It does not supersede their methods or provenance
decisions. The result remains Yachida-only and `DEVELOPMENT_ONLY`; it is not
licensed as a final abstract, Results, or publication-figure claim until Feng
and Zeller replication and the definitive three-cohort seal are complete.

## 2026-09-20 — MetaPhlAn genome-size residual audit implementation record

**Question.** Under the read-proportional reference the implanted signal is a
read fraction, so a marker-length-normalised profiler should under-report large
genomes in proportion to their length:

    log2(observed / expected_read) ~= constant - log2(G_t)

giving a slope near **-1** against log2 of the implanted target's genome size.
Replacing the implanted read fraction with the genome-equivalent implanted
signal `q_it / D_i` should remove that dependence, giving a slope near **0**.

**DECIDED: how the audit is specified.**

- New script `analysis_v2/scripts/audit_metaphlan_genome_size_residual.py`,
  wired as resumable stage 7 of `run_geff_propagation_development.sh`, reading
  only `reference_comparison/target_recovery_reference_comparison.tsv` — the
  paired primary-versus-sensitivity comparison that stage 6 already produced.
  Nothing upstream is recomputed.
- **Arms.** Primary = genome-equivalent profiler scale (expected slope 0);
  sensitivity = preserved read-proportional reference (expected slope -1). Both
  expectations are recorded as `expected_slope` and compared to the fit; neither
  is assumed, forced or rewarded.
- **Scope.** MetaPhlAn 4 only, Yachida only, `DEVELOPMENT_ONLY`. Bracken is not
  regressed, but the upstream `bracken_identical` = PASS proof is a required
  input gate read from `reference_comparison_validation.tsv`; it is never
  inferred from a filename.
- **Statistical unit: the implanted taxon.** Each target contributes one median
  `log2(observed/expected)` and the OLS is fitted across the ten target points.
  Observation rows are repeated measures of the same ten taxa; regressing them
  directly would be pseudoreplication, and no inferential standard error is
  computed from observation rows. `independent`, `community` and `pooled` are
  fitted separately; the pooled target median is taken over **all** eligible
  observations of that target, never as an average of the two population medians.
- **Uncertainty: a target-cluster percentile bootstrap.** 10,000 replicates by
  default, fixed seed **20260920**, resampling the ten targets with replacement,
  discarding any replicate with fewer than two distinct genome sizes, requiring
  at least 95% valid replicates, and taking the 2.5th/97.5th percentiles with a
  locally implemented interpolating percentile (no new SciPy dependency). A
  fresh generator with the same seed is used per (scope, arm), so the two arms
  see identical target resamples and their intervals are directly comparable.
  This is **descriptive uncertainty over the ten implanted taxa only**, stated
  in the validation output and in `DEVELOPMENT_ONLY.txt`.
- **Genome sizes.** Exclusively the FASTA-measured `target_genome_sizes.tsv`,
  revalidated against the same plausibility window `build_target_genome_sizes.py`
  enforces. No genome size is ever estimated from a recovery outcome and the
  superseded fitted 3.10 Mb constant appears nowhere
  (`fitted_genome_size_constant_used = NONE`).
- **Exact target mapping.** The comparison carries the *profiler* feature name
  while the genome-size table carries the canonical spike-panel taxon name, and
  for Fnuc those genuinely differ (`Fusobacterium nucleatum` versus
  `Fusobacterium nucleatum subsp. nucleatum`). Mapping is therefore an exact
  equality join through the frozen `examples/spike_taxon_aliases.csv` on
  `(canonical taxon name, tool)`, plus exact equality on the canonical name
  itself. The index must be injective, each label must resolve to exactly one
  observed feature, and all ten frozen labels must be present. There is no
  fuzzy, prefix, case-insensitive, normalised or substring fallback; tests
  assert that each of those is rejected.
- **Eligibility.** A ratio must be present and strictly positive and finite. No
  pseudocount is ever added. Exclusion is **paired**: an observation is dropped
  from both arms if either arm's ratio is unusable, because the arms are paired
  measurements of one physical observation and a differently filtered pair would
  make `primary_minus_sensitivity` meaningless. Exclusions are counted and
  reported per target with reasons; a target with no eligible observation in any
  scope fails the audit.
- **Fail-closed on:** missing or empty inputs; a non-empty output directory; a
  comparison lacking either profiler or either arm's columns; a blank paired
  value; a duplicate physical observation key; a missing or non-PASS Bracken
  gate; any mapping ambiguity; a missing, nonnumeric, nonfinite, non-positive,
  inconsistently duplicated or implausible genome size; a missing or unexpected
  target label; fewer than ten targets in a regression; and a bootstrap with
  fewer than 95% valid replicates.

**Status at the time of this entry:** the real Yachida result had not yet been
run. The audit was implemented and tested locally against a deterministic synthetic fixture only
(`analysis_v2/tests/test_metaphlan_genome_size_residual_audit.py`, 51 tests,
including mutation checks proving that swapping the arms or regressing
observation rows breaks named assertions). Nothing was executed on the cluster;
André runs the development pipeline personally. The audit stays
`DEVELOPMENT_ONLY` until Feng and Zeller replicate it and the three-cohort run
validates it.

**2026-09-20 addendum — four defects corrected after Codex verification.** The
estimand, target-level aggregation, paired exclusion, bootstrap, seed, scopes
and output names are unchanged. Corrected:

1. **Selected reference provenance is now carried explicitly.** The original
   row-level `reference_type` and the *selected estimand* are different things:
   a MetaPhlAn row keeps `genome_equivalent` row provenance even in the
   read-proportional arm, because that field records the source profile. The
   comparator previously carried only `reference_type`, so the sensitivity arm
   was mislabelled genome-equivalent, and the old synthetic fixture encoded that
   mistake and hid it. `compare_target_recovery_references.py` now reads
   `reference_type`, `reference_scale` and `selected_reference_type` from both
   arms' target-recovery tables and writes
   `{primary,sensitivity}_{row_reference_type,selected_reference_type,reference_scale,selected_estimand}`.
   The ambiguous `{primary,sensitivity}_reference_type` names survive only as
   backward-compatible aliases of the row type and are not read by the audit.
   Contracts, enforced in both the comparator and the audit:
   primary = `profiler_scale` / `profiler_scale_primary`; sensitivity =
   `read_proportional` / `read_proportional`; MetaPhlAn source rows must be
   `genome_equivalent` and Bracken source rows `read_proportional`; the resolved
   estimand must be `genome_equivalent` for MetaPhlAn primary and
   `read_proportional` for MetaPhlAn sensitivity and for Bracken in both arms.
   Under `profiler_scale` the selected estimand is the row's own reference,
   which is exactly why Bracken is read-proportional in both arms. In
   `metaphlan_genome_size_target_level.tsv`, `reference_type` is now the
   **selected estimand**, with `row_reference_type`, `selected_reference_type`
   and `reference_scale` kept beside it so no provenance is lost.
2. **The Yachida-only scope is proved, not asserted.** Every comparison row's
   cohort is inspected, including Bracken rows that are never regressed. The
   distinct cohort set must be exactly `{"yachida"}` — the canonical lower-case
   identifier used throughout the pipeline. Blank, Feng-only, Zeller-only and
   mixed input all fail; non-matching rows are never silently filtered.
   `observed_cohorts`, `validated_cohort_count` and `required_cohort` are
   recorded in the validation table.
3. **All physical-identity components are required and nonblank.** Every field
   of `OBSERVATION_KEY` is checked before the duplicate-key test, so two rows
   can never collide merely because a field was empty. Errors name the line and
   the field; blanks are never replaced by a sentinel; complete duplicate keys
   still fail.
4. **Absolute relative error is validated, never silently replaced.** For both
   arms of every otherwise eligible MetaPhlAn row the supplied value must be
   numeric, finite, non-negative and satisfy
   `absolute_relative_error ~= abs(observed_over_expected - 1)` within
   `abs(actual - expected) <= max(1e-12, 1e-9 * max(abs(actual), abs(expected), 1))`,
   a window sized for 17-significant-digit TSV round-tripping. Rows already
   paired-excluded for an unusable ratio never reach this check and are not
   rescued by it. Missing error fields remain structural hard failures. The
   tolerance policy is recorded in the validation metadata.

**Paired exclusion is unchanged and remains the accepted policy:** a present but
zero, negative, NaN, infinite or nonnumeric ratio in either arm excludes that
physical observation from both arms; missing or blank values are structural
failures; no pseudocounts; an entirely unusable target or scope fails; exclusion
counts and reasons are still reported.

Tests at that implementation stage: `test_metaphlan_genome_size_residual_audit.py` 80 tests and
`test_target_recovery_reference_comparison.py` 9 tests, all passing, with
mutation checks confirming each of the four gates is load-bearing. The later
21 September entry records André's completed, checksum-verified execution.

**2026-09-20 addendum 2 — stage reuse is schema-gated.** The completed Yachida
development run holds a `reference_comparison` produced before the explicit
per-arm reference-provenance columns existed. It still carries a valid
`SUCCESS`, so `run_geff_propagation_development.sh` would have reused it and
stage 7 would then have failed on the old schema. A `SUCCESS` marker is now
treated as evidence that a stage finished, not that it matches the schema the
current consumers need.

`analysis_v2/lib/stage_compatibility.sh` (`reuse_or_quarantine`) asks
`analysis_v2/scripts/check_stage_schema.py` whether a completed directory is
reusable. Headers are parsed as exact tab-separated fields inside that helper —
never an unquoted `grep` or substring test, which would accept
`primary_reference_type` where `primary_row_reference_type` is required. A
reusable `reference_comparison` needs nonempty `SUCCESS`,
`target_recovery_reference_comparison.tsv` and
`reference_comparison_validation.tsv`, the eight
`{primary,sensitivity}_{row_reference_type,selected_reference_type,reference_scale,selected_estimand}`
columns, and the `primary_selection`, `sensitivity_selection` and
`selected_estimand_contract` metrics. A reusable
`metaphlan_genome_size_residual_audit` needs its seven outputs, the
`reference_type`, `row_reference_type`, `selected_reference_type` and
`reference_scale` columns in the target-level table, and the
`selected_reference_sensitivity`, `selected_reference_primary`,
`observed_cohorts`, `absolute_relative_error_policy` and `audit_status` metrics.
An incompatible directory is moved whole to
`$RUN_ROOT/failed_attempts/<name>_<stamp>` with the reason printed, then
regenerated; both stages re-assert compatibility after regenerating. The
target-recovery arms are inputs to stage 6 and are never modified, regenerated
or deleted by this migration.

The audit additionally validates the comparator's `*_selected_estimand` columns
against the estimand derived from `reference_scale` and
`selected_reference_type`. The fields are redundant, but a contradiction means
the upstream provenance is untrustworthy, so it fails closed.

No scientific decision changed: estimand, target-level aggregation, paired
exclusion, bootstrap, seed and output names were untouched. At this point in
the chronology the execution was still pending; see the 21 September result.

## 2026-09-19 — Codex verification and exactly-one-driver release gate

**VERIFIED locally.** The numerical profiler-scale integration test, full
Python suite, DiD fixture, R parse check, and maintained shell syntax checks
pass. `collect_independent_drivers` now rejects any second direct-target row,
including an identical duplicate, so every independent physical observation
has exactly one selected-scale driver. The cluster handoff is released only for
André's `DEVELOPMENT_ONLY` execution. The continuous model has not yet run with
`mgcv`; no definitive result is licensed. The separate Figure 6 defect remains
`OPEN` below.

## 2026-09-18 — analyzer honours the selected reference scale end to end

**Problem.** `analyze_perturbation_response.py` validated profiler-scale rows
correctly but still produced **mixed-scale** predictions and superposition
results: the retained baseline was reconstructed internally as `(1-F)o`, the
response operator regressed on the implanted *read* fraction
`implanted_fraction_by_target`, and the superposition query derived the
community retained baseline as `expected - target_fraction_for_feature` while
reading read-proportional component signals.

**DECIDED and implemented.**

1. **Authoritative retained baseline.** `ResponseRow` stores the supplied
   selected retained baseline (`supplied_retained_baseline`) and
   `ResponseRow.retained_baseline` returns it verbatim. Under
   `read_proportional` that value is `(1-F)o`; under `profiler_scale` MetaPhlAn
   retains `retained_baseline_profiler_scale = (1-F)o / D_i`. The
   `(1 - effective_total_fraction) * baseline` reconstruction is removed.
   `effective_total_fraction` remains perturbation metadata and no longer
   defines any selected retained abundance.
2. **Per-perturbation operator driver.** A deterministic pre-pass
   (`collect_independent_drivers`) records, for every independent observation,
   the `direct_fraction` of its direct-implanted-target row — the implanted read
   fraction for Bracken, `q_it / D_i` for MetaPhlAn under `profiler_scale`, and
   the read fraction for both profilers under `read_proportional`. It is keyed
   by `observation_id` and cross-validated against cohort, sample, profiler and
   target. Exactly one driver per independent observation is required; missing,
   duplicate (including identical), non-positive and nonfinite drivers fail closed. That
   driver — never `row.fractions[target]` — is the regressor `x` in
   `OperatorDesign.add()`, `RegressionAggregate.add()` and the cross-cohort
   holdout prediction, which is exactly
   `selected_retained_baseline + selected_driver * fitted_slope`. Every feature
   of the perturbation, including non-implanted ones, uses that same driver.
3. **Exact superposition column selection.** `exact_superposition_rows` takes
   `reference_scale` and selects one internally consistent triple:
   `read_proportional` -> (`response_signal`, `expected_abundance_fraction`,
   `dilution_retained_baseline`); `profiler_scale` ->
   (`response_signal_profiler_scale`, `expected_abundance_profiler_scale`,
   `retained_baseline_profiler_scale`). Missing columns fail closed, so a
   profiler-scale run cannot silently fall back to the old columns. The
   retained baseline is never derived as `expected - target_fraction_for_feature`.
   The ordinal independent/community dose matching (community `dose_index - 2`)
   is unchanged.

**Unchanged by design.** Detection calculations, Bracken results under either
scale selection, and all read-proportional columns.

**Validation.** `analysis_v2/tests/test_builder_analyzer_profiler_scale.py` was
extended to a discriminating real-builder -> real-analyzer fixture: two cohorts,
both profilers, matched independent (`dose_04`) and community (`dose_06`)
perturbations, unequal target genomes (2.0 Mb / 6.0 Mb), sample-specific
`G_eff` (3.477 Mb / 4.284 Mb), and a non-implanted feature reported inside every
perturbation. Seventeen tests assert hand-computed profiler-scale numbers,
including the contents of `superposition_summary.tsv` and
`heldout_operator_validation.tsv`, that at least one MetaPhlAn result differs
between the two scale selections, that Bracken and detection are bit-identical
between them, and that corrupting the selected retained baseline or driver fails
closed. Mutation-checked: reverting any one of the three fixes above makes the
new assertions fail.

**Status.** Codex re-verification passed 19 September 2026. The cluster handoff
is released only for André's `DEVELOPMENT_ONLY` execution.
`fit_continuous_dose_response.R` was parse-checked only; `mgcv` is absent
locally and the continuous model has still never been executed.

## 2026-09-18 — PENDING (NOT fixed): Figure 6 implanted-taxon exclusion

Recorded as an open defect during the reference-scale correction and
deliberately **not** implemented in that task. Three distinct problems in the
Figure 6 community path:

- community profiles currently exclude only the **focal** target rather than all
  ten implanted taxa, so nine implanted taxa leak into the feature universe that
  is reported as non-implanted;
- community profiles may be **counted repeatedly** across `target_label` rows,
  inflating the effective number of observations;
- **feature aliases must be canonicalized before target exclusion**, otherwise an
  aliased implanted taxon survives exclusion under its source name.

Status: `OPEN`. Validation gate: a fixture in which an aliased implanted taxon
and a repeated community profile are both present must show exactly ten excluded
taxa and one counted community observation per physical profile. Until then,
Figure 6 community numbers are not licensed.

## 2026-09-04 — native profiler semantics

- Kraken2/Bracken `fraction_total_reads` and MetaPhlAn relative abundance are
  retained as profiler-native outputs; MetaPhlAn percentages are converted to
  fractions only by division by 100.
- Outputs are not claimed to measure identical cellular abundance.
- Baseline-adjusted within-profiler response is primary.
- The corrected audit passed all 720 native clean-assembly Yachida profiles.
  Definitive runs must repeat it on their exact inputs.

## 2026-09-04 — incomplete-data development run

- Ten completed strict samples were used only for integration testing while the
  prespecified 30-sample subset remained incomplete.
- The run is labelled `DEVELOPMENT_ONLY`; its estimates and p-values are not
  manuscript evidence and cannot be used for model selection.
- Its gates passed: 560 canonical rows, 500 unique native profiles, 480 paired
  positive-dose endpoints, and a converged model.

## 2026-09-04 — assembly-choice estimands

- Cleaner Pana and Pint assemblies are additive sensitivity arms; original
  results remain preserved.
- Assembly effects are estimated separately for each target and profiler.
  Equal-weight pooled effects are secondary.
- Claims concern assembly choice/quality, not a causal contamination effect,
  because strain, completeness, contiguity, and database representation differ.
- Four target-by-profiler clean-minus-original tests form the primary BH family.
  Target-specific profiler differences and pooled results are secondary.

## 2026-09-04 — repeated-measures uncertainty

- A sample random intercept alone was rejected for final assembly inference
  because it does not represent heterogeneous dose-response trajectories.
- The revised model includes sample-target-profiler random intercepts, shared
  random dose slopes, and paired clean-arm random slope deviations.
- This decision precedes the sealed 30-sample fit. Earlier random-intercept-only
  output remains diagnostic history and must not be reported.
- Literal zero p-values from numerical underflow are prohibited; extreme
  evidence is additionally represented as `-log10(p)`.

## 2026-09-04 — biological-sample-level primary assembly inference

- Development fits showed that even the random-slope GAM could give implausibly
  narrow model-based uncertainty for clean-minus-original slopes with only ten
  completed biological samples. Those development p-values are not evidence.
- The primary assembly estimator is now a two-stage paired analysis: estimate a
  six-positive-dose slope within every sample, target, profiler, and arm; then
  compare clean and original slopes within each biological sample.
- Primary uncertainty is a deterministic biological-sample bootstrap and the
  primary null test is a two-sided sample-level sign-flip test. The four
  target-by-profiler tests form one BH family.
- The GAM is retained as a secondary trajectory diagnostic. This decision was
  made before fitting the sealed 30-sample dataset and must not be changed based
  on which method produces more favorable significance.

## 2026-09-04 — paired perturbation biomarker model

- Biomarker propagation is estimated from within-sample spiked-minus-baseline
  log2 abundance changes, separately within phenotype backgrounds. Biological
  samples are the replicates; duplicated profiles are never independent.
- A fixed `1e-8` fraction pseudocount replaces outcome-dependent minimum-value
  pseudocounts. The off-target feature universe is fixed across all doses and
  requires 10% nonzero prevalence; the intended target is always retained.
- BH correction is performed across species within each target, arm, profiler,
  background, and dose. Positive calls at q <= 0.05 are primary; q <= 0.10 is a
  prespecified sensitivity analysis.
- This controlled perturbation contrast is not called CRC-versus-control.
  Baseline disease contrasts require a separate, still-unfrozen model.

## 2026-09-04 — native disease-biomarker model

- Actual disease contrasts are fitted separately by cohort, target, assembly
  arm, profiler, and dose; cohorts are not pooled for primary inference.
- CRC versus Control is primary and Adenoma versus Control is secondary. The
  primary adjustment set is age and sex; adding BMI on complete cases is a
  prespecified sensitivity analysis.
- Native abundance fractions are transformed as `log2(x + 1e-8)`. HC3 robust
  standard errors are used. The cross-dose species universe requires 10%
  prevalence, with the intended target always retained.
- BH correction is across species within the exact cohort, population, target,
  assembly, profiler, dose, contrast, and model context. Baseline disease calls
  are observed calls, so baseline-to-dose Jaccard stability is meaningful.
- Development data may exercise this model but cannot support manuscript
  estimates. Definitive claims require sealed complete cohort inputs.
- Covariates with no variation in an analysis population are non-estimable and
  cannot confound within that population. They are omitted automatically and
  recorded per result; other rank deficiency is not silently repaired.

## 2026-09-05 — cross-cohort synthesis

- Cohort-specific estimates remain primary; no sample-level pooled analysis
  that ignores cohort is permitted.
- Shared disease-biomarker effects use REML random effects with conservative
  modified Hartung–Knapp uncertainty. Heterogeneity, prediction intervals,
  direction agreement, and leave-one-cohort-out estimates are mandatory.
- Only features with positive finite uncertainty in all three prespecified
  cohorts are pooled. Missing coverage remains explicit.
- With three cohorts, heterogeneity is imprecise; synthesis supports cautious
  generalization rather than proof of universal consistency.

## 2026-09-05 — disease-biomarker propagation semantics

- An organism implanted across phenotype groups is not expected to become a
  disease biomarker. Its disease significance is a spurious-association
  diagnostic, not target recall.
- **Superseded 2026-09-14:** the target can already carry a baseline disease
  contrast. Post-spike significance alone is therefore neither recall nor a
  spurious-association label; use baseline-to-dose effect change and the
  ideal-reference residual instead.
- Disease propagation uses two-sided significant call sets and reports retained,
  lost, and gained biomarkers, baseline retention, Jaccard stability, direction
  flips, and changes in baseline-biomarker effects.
- Empty-versus-empty biomarker sets have undefined Jaccard similarity. They are
  not assigned a value of one, which would imply evidence of perfect stability.

## 2026-09-08 — frozen baseline disease feature universe

- Disease-biomarker eligibility is determined once from unmodified baseline
  profiles within each cohort and profiler, using 10% nonzero prevalence.
- All prespecified implanted targets are added together as a common diagnostic
  set. A target or assembly arm cannot change the multiplicity burden.
- Baseline coefficients and BH-adjusted q-values are computed once and reused
  across every target, arm, and dose. Positive-dose models use the identical
  feature universe, preventing post-perturbation feature selection.

## 2026-09-08 — direct artificial-biomarker recovery report

- Spiked versus matched unmodified libraries form a separate, direct recovery
  analysis in which the implanted target is a prespecified artificial biomarker.
- Primary outputs are target recall, precision, off-target enriched calls, the
  paired target effect, and the minimum tested dose producing target recovery at
  BH q <= 0.05. The q <= 0.10 result is a labelled sensitivity analysis.
- Minimum detected dose is restricted to the tested dose grid and is not called
  a continuous detection limit. Results remain separate by phenotype background.
- Exact achieved fractions are retained as evidence, but aggregation uses the
  six frozen nominal dose levels so integer read allocation cannot fragment a
  single experimental dose into several summary groups.
- Reports validate population-specific target-dose grids: six target fractions
  for independent spikes and seven target fractions for the ten-member legacy
  community mixture. Total community dose is never mistaken for per-target dose.
- This analysis complements rather than replaces disease-versus-control
  propagation, and its effects describe sequencing evidence rather than cells.

## 2026-09-08 — calibration-to-biomarker linkage

- Sealed quantitative response endpoints are joined to sealed artificial-
  biomarker calls; neither source analysis is refitted by the linkage module.
- Continuous read-perturbation response ratio is primary. A median ratio below
  0.8, within 0.8--1.2, or above 1.2 is labelled under-response,
  read-proportional-band, or over-response only for descriptive presentation.
- Associations with target effects, recall, precision, and off-target burden
  are descriptive and cannot establish that calibration error causes a
  biomarker result. Profiler-native response is not cellular-abundance accuracy.

## 2026-09-08 — pooled primary artificial-biomarker inference

- The primary direct spiked-versus-baseline analysis pools all eligible paired
  biological samples across phenotype backgrounds; it is not a disease contrast.
- Control-, Adenoma-, and CRC-stratified fits remain secondary.

## Pooled artificial-biomarker estimand and off-target identities

- The artificial biomarker created by controlled read implantation is tested
  primarily by pooling biological samples and pairing each spiked profile with
  its own unmodified profile.
- Phenotype-background analyses are secondary heterogeneity/descriptive
  analyses, not substitutes for the pooled primary estimand.
- Detection thresholds include first observed and first sustained detection.
- Off-target discovery is recorded both as a burden and as a taxon-level
  recurrence ledger so cross-mapping or compositional effects can be audited.
- Calibration linkage uses profiler-native baseline-adjusted changes and keeps
  pooled and phenotype-stratified estimands explicitly separated.

## Definitive Yachida execution gate

- Manuscript-labelled Yachida results require a checksummed 201-sample
  production seal and a validated canonical table covering all community
  samples and the nested 30-sample independent subset.
- The definitive driver is fail-closed and cannot inherit development status.
- Independent and community calibration models are both mandatory. The pooled
  paired artificial-biomarker analysis is primary; phenotype-stratified and
  clean-assembly analyses remain secondary.
- Community target doses are reconstructed from the exact integer allocation
  used during spike generation. The target-specific read fraction is its
  allocated pair count divided by the final library pair count; equal division
  of the nominal community fraction is not assumed.

## Completed assembly-choice sensitivity endpoints

- Pana/Pint clean-versus-original detection uses paired binary outcomes within
  sample, target, profiler, and dose; exact McNemar p-values are BH-adjusted
  across the 24 secondary comparisons.
- Pooled paired artificial-biomarker models are applied identically to both
  assembly arms. Target calls/effects, precision, and off-target burden are
  compared descriptively, with complete off-target taxon ledgers.
- These results are labelled assembly-choice sensitivity and cannot isolate a
  causal effect of contamination from strain or database representation.

## Shared CRC-cohort definitive framework

- Feng and Zeller share implementation code but remain separate studies for
  estimation and reporting. Their production populations contain 154 and 156
  eligible samples, respectively; each individual-spike subset is the frozen
  balanced 10/10/10 selection.
- A cohort cannot enter definitive analysis without a checksummed upstream
  seal, complete sample-flow ledger, covariate-missingness audit, exact expected
  profile counts, validated canonical table, and native-profile audit.
- Three-cohort synthesis is allowed only from definitive Yachida, Feng, and
  Zeller packages; incomplete feature coverage remains in the synthesis ledger.

## 2026-09-09 — Cross-cutting analysis policy frozen

- Native abundance greater than zero is the operational detection definition;
  it is not interpreted as a common analytical detection limit or biological
  presence/absence truth.
- Unconditional quantitative summaries retain non-detections as zero. Paired
  calibration endpoints use no pseudocount or transform, retain negative
  baseline-adjusted responses, and define response ratios only for positive
  implanted fractions.
- The fixed `1e-8` fraction pseudocount applies only to the prespecified log2
  biomarker models. Detected-only and species-closed results are secondary.
- `ANALYSIS_POLICY.tsv` records these rules and all BH-family boundaries.
  Definitive runners validate and checksum this policy into every cohort
  package before analysis begins.

## 2026-09-14 — artefact filtering and biomarker-fragility interpretation

- Loss of a native disease biomarker after controlled implantation is defined
  as **perturbation fragility**, not evidence that the baseline association is
  a false positive.
- The implanted target and its frozen profiler-specific alias are excluded from
  the primary **bystander-distortion** analysis. Direct-target recovery and
  attenuation remain separate mechanistic endpoints.
- Community and independent baselines are population-specific unspiked panels;
  differences in their baseline call counts are not effects of implantation.
- A development-only, blind static-blacklist screen transferred artefact scores
  between Feng and Zeller in both directions without protecting held-out
  targets. Most candidate rules removed many off-target calls only by also
  losing substantial target recall. No filtering threshold was selected.
- No manuscript claim that spike-derived filtering removes false positives is
  permitted until a rule is frozen in training data, applied before model
  fitting, followed by complete model and BH-family refitting, and evaluated in
  held-out data with both off-target reduction and target-recall loss reported.
- Native disease-effect distortion will use an ideal read-proportional
  counterfactual: `e = (1 - F)o` for bystanders and `e = (1 - F)o + f` for the
  implanted target. The primary paired difference-in-differences estimand is
  the phenotype effect on observed-versus-ideal residual abundance. Continuous
  effect distortion is primary; significance-state transitions are secondary.
- The hypothesis that fragile biomarkers are enriched for false positives will
  be tested against null-label discoveries and externally directionally
  replicated biomarkers. Fragility becomes a filter only if discrimination and
  recall preservation succeed out of cohort; otherwise it remains a robustness
  annotation.

## 2026-09-14 — revised analysis policy version 2

- Policy version 2 supersedes version 1 only for the ongoing revised analysis;
  it does not alter the tagged bioRxiv-v1 workflow.
- The primary native robustness scope is explicitly limited to non-target
  bystanders. Direct-target behavior is a secondary mechanistic endpoint.
- The ideal read-proportional reference model treats continuous distortion as
  primary and adjusts across the frozen feature universe within each cohort,
  study population, target, assembly arm, profiler, and dose.
- No global claim that any implantation causes distortion is licensed by the
  context-wise families. Such a claim requires a separately specified
  hierarchical multiplicity procedure.

## 2026-09-14 — cross-cohort abundance-calibration experiment

- Taxonomic distortion is evaluated before biomarker selection using paired
  observed-versus-ideal abundance residuals across the complete frozen
  baseline-prevalence universe.
- Calibration coefficients are learned without the validation cohort. Only
  positive residual slopes supported by at least 10 biological samples and a
  60% positive-residual fraction are eligible for conservative subtraction.
- Independent experiments retain target-specific coefficients. Community
  experiments use one joint-mixture coefficient because component fractions
  are not independently identifiable in the fixed mixture design.
- Direct implanted target features are never subtracted. Corrected abundances
  are bounded below by zero; unsupported and negative corrections are zero.
- Validation requires complete refitting of the unchanged disease model and
  BH families. Benefits (abundance/effect restoration, rescued baseline calls,
  removed induced calls) and harms (damaged baseline calls and created calls)
  are reported together. No automatic feature removal is authorized.

## 2026-09-14 — achieved-dose reporting tolerance

- Reporting maps achieved read fractions to the nearest population-specific
  frozen nominal dose with a 5% relative tolerance. Integer implanted-read
  counts produce several-percent relative deviations at the smallest community
  doses; the former 0.1% tolerance incorrectly rejected valid observations.
- Five percent remains far below the separation between adjacent frozen doses,
  so every accepted value has an unambiguous nominal level. Exact achieved
  fractions remain in the detailed evidence tables.

## 2026-09-17 — profiler-specific expected response for MetaPhlAn

- **Status: DECIDED IN PRINCIPLE; implementation and `G_eff` validation are
  pending.** The earlier use of `e = (1-F)o + f` for both profilers is retained
  only as an explicitly named read-proportional sensitivity analysis.
- Bracken `fraction_total_reads` retains the read-proportional primary
  reference.
- MetaPhlAn `relative_abundance` is a marker-length-normalized,
  genome-equivalent-like composition. Paired baseline subtraction does not by
  itself convert an implanted read fraction into this estimand.
- For target assembly length `G_t` and an independently derived effective
  baseline community genome size `G_eff,i`, define
  `q_it = f_it * G_eff,i / G_t` and
  `Q_i = sum_k(f_ik * G_eff,i / G_k)` across implanted members. The first-order
  MetaPhlAn expectation is
  `((1-F_i)o_it + q_it) / ((1-F_i)+Q_i)` for an implanted target and
  `(1-F_i)o_ij / ((1-F_i)+Q_i)` for a non-implanted feature.
- `G_t` must come from the exact checksummed FASTA used for spike simulation.
  `G_eff,i` must be frozen independently of the observed post-spike recovery;
  choosing it to make MetaPhlAn appear unbiased is prohibited. A pooled fitted
  value is diagnostic only.
- A mapped-community fraction correction (`f/m`) is not interchangeable with
  genome-equivalent scaling and is not the selected solution.
- MetaPhlAn quantitative recovery, response residuals, recovery classes,
  cross-talk residuals, community superposition, and ideal-counterfactual
  distortion remain provisional until the corrected fields and tests are
  implemented. Detection, native observed profiles, baseline disease models,
  and baseline cross-cohort biomarker replication do not require re-profiling.
- The correction is downstream-only using existing exact fractions, FASTAs,
  and retained profiler outputs. Raw MetaPhlAn mapouts were normally temporary
  because production used `KEEP_MAPOUT=0`; an estimator-exact marker-level
  reconstruction is therefore optional and would require remapping.
- Development evidence linking MetaPhlAn read-reference error to target genome
  size is hypothesis-generating until reproduced at biological-sample and
  cohort level with authoritative assembly lengths and uncertainty.

## 2026-09-17 — G_eff source and profiler-scale endpoint implementation

- **DECIDED.** Per-SGB genome sizes `G_j` for the effective community genome
  size `G_eff,i` are taken from the **MetaPhlAn vJan25 database itself**. It is
  the reference system MetaPhlAn normalizes against, so the correction stays
  internally consistent and needs no lossy SGB-to-species remapping. UHGG
  representative lengths, NCBI per-species lengths, and a cohort-level constant
  were considered and rejected as primary; any of them may be reported as a
  prespecified sensitivity. The extraction from the database and the
  mapping-coverage audit on sealed Yachida baselines remain outstanding.
- **DECIDED.** New profiler-scale fields are **additive**. Existing
  read-proportional field names are unchanged, so no downstream consumer breaks
  and migration is incremental. The earlier proposal to rename them to
  `*_read_reference` is not adopted; `reference_type` carries the distinction
  explicitly on every row.
- **IMPLEMENTED** in `scripts/derive_paired_endpoints.py`. Bracken keeps
  `e = (1-F)o + f`. MetaPhlAn uses `q_it = f_it * G_eff,i / G_t`,
  `Q_i = sum_k q_ik`, `e = ((1-F)o + q)/((1-F) + Q)`, with matching retained
  baseline, implanted signal, recovered signal, response ratio, and residual.
  Eleven fields were added, led by `reference_type`.
- `Q_i` is accumulated per `profile_id`, which separates the designs without a
  special case: an independent spike contributes one member and a community
  spike contributes all ten.
- Fail-closed gates: MetaPhlAn rows present without an explicit
  `--metaphlan-reference` choice; `genome_equivalent` without both genome
  tables; any unmapped target or sample; and community member fractions that
  do not reproduce `spike_fraction_total` within 1% relative, which would mean
  the canonical table is missing members and `Q_i` is understated. There is no
  default reference scale, because defaulting is how the read-fraction versus
  genome-equivalent estimand mismatch entered the analysis.
- Fixture tests in `tests/test_paired_endpoints_profiler_scale.py` cover all
  four gates, equal and unequal genome sizes, community summation, zero
  baseline, Bracken invariance, and the `read_proportional` opt-out.
  `test_canonical_input.py` and `test_paired_endpoints.py` still pass.
- Still outstanding: the `G_eff,i` derivation itself, the non-implanted
  reference `e_MP,ij` in the atlas input, and the remaining consumers
  (`build_perturbation_response_input.py`, `analyze_perturbation_response.py`,
  `summarize_target_recovery.py`, `fit_continuous_dose_response.R`,
  `fit_ideal_counterfactual_did.R`). No definitive MetaPhlAn recovery claim is
  licensed until those are complete and validated on sealed Yachida.

## 2026-09-17 — effective community genome size G_eff,i (PROVISIONAL)

- **PROVISIONAL, not DECIDED.** The production MetaPhlAn database was not
  available in the implementation environment, so the authoritative genome-size
  mapping was not extracted and the sealed Yachida coverage audit was not run.
  The policy below is a stated starting point. It becomes `DECIDED` only after
  the cluster steps in `METAPHLAN_EFFECTIVE_GENOME_SIZE.md` succeed and the
  observed coverage distribution is inspected.
- **Estimand.** Sample-specific `G_eff,i = sum_j(c_ij * G_j) / sum_j(c_ij)` over
  eligible mapped species-rank features of the unspiked baseline profile.
  Sample-specific is primary; cohort-median and global constants are labelled
  sensitivities only.
- **`G_j` source** is the production MetaPhlAn vJan25 database taxonomy mapping,
  per the 2026-09-17 decision above. Marker length is never a size source; a
  test asserts a marker record cannot enter the output.
- **Eligibility.** Species-rank (`s__`) terminal clades only. Higher ranks are
  excluded because MetaPhlAn repeats mass at every rank; `t__` rows are excluded
  because they are children of species rows and would double-count. Excluded and
  unclassified mass is reported separately, never silently folded in.
- **Unmapped features are never imputed.** They are dropped from numerator and
  denominator and reported, and `mapping_coverage` states how much of the
  community the estimate rests on.
- **Provisional thresholds.** 95% minimum mapped eligible abundance as primary,
  90% as a prespecified sensitivity, zero tolerated exclusions by default.
  These must be re-examined against real Yachida coverage before freezing, but
  must never be chosen by comparing recovery performance.
- **Independence from spike outcomes is enforced, not merely asserted.** The
  manifest must declare `spike_fraction_total == 0` and
  `profile_id == baseline_profile_id` for every row, using the canonical
  contract's own invariants rather than a filename heuristic; there is no
  override flag. A test writes a spiked profile with very different composition
  on disk and asserts `G_eff,i` is unchanged.
- The provisional 3.10 Mb value fitted from development recovery data is **not
  used**. It was derived from post-spike outcomes, which this policy forbids.
- **Implemented:** `scripts/extract_metaphlan_genome_sizes.py` (schema-tolerant
  database extraction with provenance and checksums; fails closed on an
  unrecognised schema and prints the observed schema) and
  `scripts/compute_effective_genome_size.py` (per-profile `G_eff`, coverage,
  exclusions, unmapped-feature report, audit, checksums, deterministic sorted
  output). Documented in `METAPHLAN_EFFECTIVE_GENOME_SIZE.md`.
- **Tested:** `tests/test_effective_genome_size.py` (13 cases) and
  `tests/test_metaphlan_genome_size_extraction.py` (6 cases).
- **KNOWN INCOMPATIBILITY, not yet fixed.** `derive_paired_endpoints.py` now
  requires `--metaphlan-reference` whenever MetaPhlAn rows are present. These
  four runners still invoke it without that argument and will fail closed on any
  MetaPhlAn-containing input: `run_yachida_definitive_analysis.sh`,
  `run_crc_cohort_definitive_analysis.sh`, `run_assembly_sensitivity.sh`, and
  `prepare_three_cohort_development_input.sh`. This is deliberate — failing is
  correct until each runner states its reference scale — but it must be resolved
  before any definitive run.

## 2026-09-17 (corrective) — G_eff rank, input enforcement, provenance (still PROVISIONAL)

Corrects the entry above after independent review. Status remains
**PROVISIONAL**; nothing here has been run against the real database.

- **Rank corrected — this was a genuine defect.** The first implementation
  weighted `s__` rows and excluded all `t__` rows. The production vJan25
  taxonomy is SGB-based, with keys ending `...|s__Species|t__SGB123`, evidenced
  by `scripts/build_reference_representation_table.py::metaphlan_sgb_map()`
  requiring both `s__` and `t__SGB\d+` in one key, and by
  `workflows/metaphlan4/postprocess_local.sh` existing to strip `t__` rows from
  raw profiles. `source_profile` is the raw `<sample>.metaphlan.tsv`, so SGB
  rows are present. Against the real database the original code would have
  mapped nothing and reported near-zero coverage. Earlier fixtures hid this by
  using species-terminal synthetic keys.
- **`--profile-rank sgb` is now the default:** terminal rows at or below species
  rank, keyed on full lineage for exact database matching. Terminal-only
  eligibility is what prevents double counting, since a species row with an SGB
  child is non-terminal. A test asserts eligible mass is 100, not 200, when
  species and SGB rows each carry the same mass. A terminal `s__` row with no
  SGB child stays eligible and is reported unmapped, lowering coverage honestly.
- `--profile-rank species` is retained as a documented alternative. No
  SGB-to-species aggregation policy is implemented, and none is invented
  silently.
- **Rank incompatibility fails closed** in both directions, tested. This is what
  prevents a silent zero-coverage run.
- **MetaPhlAn-only enforcement** is now inside the script: `profiler` is a
  required column and every included row must be `metaphlan4`. Mixed and
  Bracken-only manifests fail. Previously a Bracken table could have been parsed
  as a MetaPhlAn profile.
- **Inclusion policy:** `include`/`exclusion_reason` required. Default
  `require_included_only` fails if any `include == 0` row is present;
  `filter_excluded` drops them with a complete audit. `include == 0` without a
  reason always fails.
- **Repeated baseline rows** are collapsed only after verifying agreement on
  nine invariant fields; conflicts fail and name the field. Previously
  `setdefault` silently kept the first row.
- **Source-profile provenance:** `source_profile_checksums.tsv` records cohort,
  sample, profile, resolved path, size, and SHA-256, hashing each physical file
  once after dedup validation, and is included in the checksum manifest.
- **Runbook repaired:** positional `awk -F '\t' '$13 == 0'` is replaced by
  `scripts/select_baseline_manifest.py`, which selects by column name, enforces
  profiler/dose/identity/inclusion, validates repeated rows, and writes its own
  audit and checksums. The runbook now uses exact mapping paths rather than a
  wildcard, marks outputs `DEVELOPMENT_ONLY`, and separates the 95% primary from
  the 90% sensitivity.
- **Tests:** `test_effective_genome_size.py` rewritten with hierarchical
  SGB-terminal fixtures (17 cases), `test_select_baseline_manifest.py` added (6
  cases), extractor fixtures corrected to SGB-terminal keys.
- **Documentation claims were reduced** to what the tests establish.
- **Unchanged and still outstanding:** the real database has not been inspected,
  the sealed Yachida coverage audit has not run, and the true SGB mapping
  coverage is unknown. The known runner incompatibility with
  `derive_paired_endpoints.py` also remains unresolved.

## 2026-09-17 (corrective 2) — G_eff is sample-wide, not per-population (PROVISIONAL)

Resolves the canonical population semantics completely. Status remains
**PROVISIONAL**: no cluster access was available, so the real database and
sealed Yachida audit still have not been run.

- **Canonical fact.** `build_crc_cohort_canonical_input.py` writes the *same
  physical* zero-dose MetaPhlAn profile for both analytical populations. Both
  baseline rows come from
  `make_row(meta, population, target, profiler, sample, sample, 0.0, 0.0, 0, ...)`
  against the same `unique_file(baseline_root, sample + ".metaphlan.tsv")`, so
  they share `profile_id == sample_id`, `baseline_profile_id == sample_id` and
  one `source_profile`, and differ only in `analysis_population`.
- **DECIDED (implementation semantics, not the scientific policy).** `G_eff,i`
  is a property of that single physical baseline. Identity is
  `(cohort, sample_id)`, one row is emitted per biological sample, and
  `analysis_population` is written **blank**. The blank value is functional, not
  cosmetic: `derive_paired_endpoints.py` resolves
  `(cohort, sample_id, analysis_population)` and falls back to
  `(cohort, sample_id, "")`, so one sample-wide value serves both populations.
  Neither population label is arbitrarily retained.
- **Defect fixed.** The previous version keyed dedup on bare `profile_id` and
  treated `analysis_population` as an invariant, so real canonical input would
  have failed with a spurious conflict. `profile_id` also equals the sample id,
  so it would have collided across cohorts. Keys now always include cohort.
- **Selector** now filters on cohort, `profiler == metaphlan4`, both spike
  fractions zero, `profile_id == baseline_profile_id` and `include == 1`;
  collapses on `(cohort, sample_id)`; blanks the population; writes a per-row
  exclusion ledger; requires a reason on every dropped `include == 0` row; and
  enforces `--expected-profiles N` by withholding `SUCCESS` and exiting non-zero.
- **Checksum ownership** is validated: one resolved path may be hashed once, but
  a path claimed by two different `(cohort, sample_id)` identities now fails
  instead of being absorbed by `setdefault`.
- **Integration test added** — `tests/test_geff_canonical_integration.py` builds
  a fixture mirroring the producer (two samples, one in both populations, both
  profilers, per-target baseline repetition, positive doses) and runs canonical
  -> selector -> G_eff -> `derive_paired_endpoints --metaphlan-reference
  genome_equivalent`. It asserts one G_eff row per sample, blank population,
  both populations resolving the same `G_eff` and the same `q_it`, expected-count
  enforcement in both directions, conflicting-baseline failure, cross-cohort
  distinctness, and no species/SGB double counting. Its schema is imported from
  the canonical test helper so it cannot drift from the producer contract.
- **Runbook** now passes `--cohort yachida --profiler metaphlan4
  --expected-profiles 201` and independently verifies the row count with awk and
  the presence of `t__SGB` rows in a real profile, rather than trusting prose.
- **Still outstanding:** vJan25 schema inspection, authoritative mapping
  extraction, the 201-profile selection against sealed Yachida, real SGB
  coverage, and the 95%/90% audits. The `derive_paired_endpoints.py` runner
  incompatibility also remains unresolved.

## 2026-09-17 — G_eff derivation validated on sealed Yachida (DECIDED)

This entry supersedes the provisional status in the preceding G_eff entries.
André ran the cluster audit manually; no agent accessed or operated the cluster.

- **Production database:** MetaPhlAn 4.2.2 vJan25 pickle
  `mpa_vJan25_CHOCOPhlAnSGB_202503.pkl`, SHA-256
  `e7d23a73a7959b4f41af0bbe403f4b5bbb7c1879528d376d146ee5294515df9a`,
  163,915,626 bytes. Inspection confirmed taxonomy keys ending in
  `t__SGB...` and tuple values containing plausible whole-genome lengths;
  marker lengths are stored separately.
- **Authoritative mapping:** 58,331 taxonomy entries yielded 58,216 usable
  terminal-SGB genome lengths; 115 entries were skipped. All extracted ranks
  were `t`; lengths ranged from 121,861 to 49,705,299 bp.
- **Sealed input:** canonical Yachida input SHA-256
  `251d0ed2df12d49cad12907d66e1273e29c9416de5a2223ca613756b0c4434cd`,
  36,360 included rows and 201 samples. The selector produced exactly one
  sample-wide baseline for each of 201 samples; 30 samples occurred in both
  analytical populations and correctly collapsed to one physical baseline.
- **Coverage:** all 201 raw baseline profiles contained SGB rows. Exact-lineage
  mapped-abundance coverage had min/median/max 1.0, with zero unmapped features
  and zero excluded samples. The prespecified 95% primary and 90% sensitivity
  analyses therefore both passed and produced identical included sets.
- **Observed reference distribution:** sample-specific `G_eff` ranged
  2,630,057.998--4,844,354.509 bp, median 3,477,488.241 bp.
- **Frozen policy:** one sample-wide `G_eff` per `(cohort, sample_id)`, derived
  exclusively from the unspiked native MetaPhlAn profile as the abundance-
  weighted mean production-database genome length over exactly matched
  terminal-SGB lineages. Unmapped mass is not imputed. Primary minimum coverage
  remains 95%, with 90% as the prespecified sensitivity.
- **Audit bundle:**
  `work/yachida_geff_audit_20260917T144215Z`.
- **Scope boundary:** this decision freezes only the derivation of `G_eff`.
  Downstream consumers and definitive runners still require migration and
  sealed-Yachida validation. No final MetaPhlAn recovery, cross-talk,
  superposition, or distortion claim is licensed yet.

## 2026-09-17 (propagation, partial) — runners migrated; consumers BLOCKED

The G_eff policy is **DECIDED** (sealed-Yachida audit run manually by André:
201/201 baselines, 58,216 terminal-SGB lengths, coverage 1.0 min/median/max,
zero unmapped, zero excluded, both gates passed, G_eff 2.630-4.844 Mb, median
3.477 Mb, bundle `work/yachida_geff_audit_20260917T144215Z`, database SHA-256
`e7d23a73a7959b4f41af0bbe403f4b5bbb7c1879528d376d146ee5294515df9a`). Nothing in
this entry changes it. The previously fitted 3.10 Mb value remains unused.

**Implemented — runner migration (requirement B).**
- New `analysis_v2/lib/metaphlan_reference.sh`, sourced by all four runners:
  `run_yachida_definitive_analysis.sh`, `run_crc_cohort_definitive_analysis.sh`,
  `run_assembly_sensitivity.sh`, `prepare_three_cohort_development_input.sh`.
- `derive_endpoints_with_references` detects MetaPhlAn rows by **column name**
  in the canonical input, then requires `TARGET_GENOME_SIZES` and
  `EFFECTIVE_GENOME_SIZES`, verifies both exist and are non-empty, and fails
  closed otherwise. There is no default reference scale.
- Primary genome-equivalent output goes to `<root>/endpoints`; the labelled
  read-proportional sensitivity goes to
  `<root>/endpoints_read_reference_sensitivity`, marked `SENSITIVITY_ONLY.txt`.
  The helper refuses to run if either directory already exists, so a
  sensitivity result can never overwrite or be confused with the primary.
- Both arms write `metaphlan_reference_provenance.tsv` recording role, resolved
  paths, input SHA-256s, canonical input and UTC timestamp. No dated Yachida
  audit directory is hard-coded; paths are supplied by environment variable and
  their resolved values are recorded.
- Bracken-only canonical inputs take the original code path unchanged.
- Tested by `tests/test_metaphlan_reference_runner.py`.

**BLOCKED — non-implanted reference (C) and downstream consumers (D).**
Two stop conditions in the assignment are met, so this was reported rather than
improvised:

1. *Uncommitted work overlaps incompatibly.* Five of the consumers to migrate
   are **untracked** in-flight files with no git history:
   `build_perturbation_response_input.py` (531 lines),
   `analyze_perturbation_response.py` (1,481),
   `summarize_target_recovery.py` (222),
   `fit_ideal_counterfactual_did.R` (522),
   `build_perturbation_reliability_scores.py`. Rewriting the expected-abundance
   core of unversioned work offers no recovery path if it is wrong.
2. *A consumer is already failing.* `tests/test_target_recovery_summary.py`
   fails **at baseline**, before any change in this session
   (`detected_samples` 0 != 1). Migrating a consumer whose own test is red would
   confound that failure with the reference change.

A third, practical constraint: `fit_continuous_dose_response.R` and
`fit_ideal_counterfactual_did.R` need `mgcv`/`sandwich`, and Parquet reading
needs `arrow`; none are installed locally (they live in the frozen analysis
image), so their reference-type behaviour cannot be regression-tested here.

Unmodified originals of all five untracked consumers were copied to the session
scratchpad before any work began; none were edited.

**What C will require when unblocked.** `build_perturbation_response_input.py`
already carries `implanted_fraction_by_target` and
`effective_total_fraction` per observation, so `Q_i` is reconstructible without
schema change. The work is to join `G_eff` on `(cohort, sample_id)` and `G_t` on
target label, then add genome-equivalent variants of
`dilution_retained_baseline`, `expected_abundance_fraction`, `response_signal`,
`response_delta`, the bounded errors and `quantitative_log2_ratio`, plus a
`reference_type` column — additively, never overwriting the read-reference
fields, whose semantics downstream code depends on.

## 2026-09-17 (propagation in progress — SUPERSEDED, see the corrective entry below) — genome-equivalent reference through the pipeline

Supersedes the "propagation, partial" entry above. The G_eff policy is unchanged
and remains DECIDED; the fitted 3.10 Mb value is still unused.

**Invalid blocker corrected.** The previous entry treated untracked consumers as
a blocker. That was wrong: untracked work can be safely extended once it is
backed up and audited. All modified untracked files were SHA-256'd and copied
byte-for-byte outside the repository before editing, with an audit recording
path, backup path, hash, size and mtime. Diffs against those originals are
provided. The pre-existing red test was likewise not a blocker.

**Baseline test repaired.** `test_target_recovery_summary.py` failed because its
fixture contained only `analysis_population='independent'` while
`fnuc_baseline_prevalence.tsv` deliberately restricts to `community`, so
`detected_samples` summed to 0. The **fixture** was repaired by adding two
community rows (one detected at baseline, one not); the production population
filter was left untouched because it is correct. The test now also asserts the
prevalence output is community-only.

**Non-implanted MetaPhlAn reference implemented** in
`build_perturbation_response_input.py`:
- `Q_i = sum_k(f_ik * G_eff,i / G_k)` is computed once per perturbed profile in
  `observation_genome_equivalents`, over **every** implanted member;
- implanted target: `e = ((1-F)o + q)/((1-F) + Q)`;
- non-implanted feature: `e = ((1-F)o)/((1-F) + Q)`;
- Bracken is untouched and keeps `(1-F)o + f`.
Fails closed on: MetaPhlAn rows without an explicit `--metaphlan-reference`;
missing genome-size or G_eff tables; an implanted target lacking `G_k`; a
MetaPhlAn sample lacking `G_eff`; incomplete community membership; a
non-positive renormalisation denominator; and — newly enforced — implanted
member fractions that do not reconstruct the recorded total within the
prespecified 5% relative tolerance. That last check previously existed only as a
reported number, never a gate.

**Field semantics.** `expected_abundance_fraction` keeps its original
read-proportional meaning for every profiler and is now accompanied by the
explicit alias `read_proportional_reference`. The primary scale lives in new
columns: `reference_type`, `expected_abundance_profiler_scale`,
`retained_baseline_profiler_scale`, `response_signal_profiler_scale`,
`response_delta_profiler_scale`, `quantitative_log2_ratio_profiler_scale`,
`effective_community_genome_size_bp`, `target_genome_size_bp`,
`implanted_genome_equivalent_fraction`,
`total_implanted_genome_equivalent_fraction`. No existing field changed meaning.

**Consumers migrated.** New `scripts/reference_scale.py` centralises selection.
`summarize_target_recovery.py` and `analyze_perturbation_response.py` now take a
**required** `--reference-scale {profiler_scale|read_proportional}` and alias the
chosen reference onto the canonical column names, so their downstream logic and
Bracken numbers are unchanged. Recovery outputs carry `profiler`,
`reference_scale` and `reference_type`. Mixed reference types are rejected
unless the caller groups by `reference_type`, and requesting `profiler_scale`
against a table lacking those columns fails loudly (a defect found by the new
tests: the TSV read path originally bypassed the remap and fell back silently).

**R consumers.** `fit_ideal_counterfactual_did.R` and
`fit_continuous_dose_response.R` now require an explicit `--reference-scale`.
`read_proportional` is numerically unchanged; `profiler_scale` **fails closed**
because the DiD reconstructs its ideal reference from the focal target alone and
therefore cannot form `Q_i` for community perturbations. Implementing it there
requires a per-member implanted-fraction input that its manifest does not carry.
This is recorded as the one remaining consumer gap rather than approximated.

**Runner helper hardened.** A canonical table without a `profiler` column now
fails instead of being treated as Bracken-only; reference tables are validated by
schema, not just non-emptiness; every MetaPhlAn cohort/sample must appear in the
G_eff table; the canonical checksum is recorded; a missing `SUCCESS` in either
arm aborts so a partial run cannot look complete; and disabling the sensitivity
is recorded in the primary provenance.

**Tests.** Twelve Python suites and the R DiD suite pass, including new
`test_perturbation_response_genome_equivalent.py` (11 cases) and an end-to-end
chain in `test_geff_canonical_integration.py`: canonical -> baseline selection ->
G_eff -> paired endpoints -> response input -> target-recovery summary.
`mgcv`, `sandwich` and `arrow` are absent locally and no image is built, so
`fit_continuous_dose_response.R` was validated by R `parse()` only.


## 2026-09-17 (corrective) — profiler-aware validation, implanted signal, R model paths

Supersedes the entry above, which claimed completeness prematurely. The G_eff
policy is untouched and remains DECIDED.

**Defect fixed — the validator rejected valid primary tables.** A primary table
normally and correctly holds two row-level reference types at once
(`kraken2_bracken -> read_proportional`, `metaphlan4 -> genome_equivalent`). The
previous `require_single_reference_type` rejected exactly that. Validation is now
**within profiler and analysis arm**, never globally: a profiler carrying more
than one reference type fails, Bracken must be `read_proportional`, MetaPhlAn
must be `genome_equivalent`, unknown profilers and blank values fail.
`read_proportional` selects the preserved read fields from the same primary
table and reports `selected_reference_type = read_proportional`; row-level
provenance is never rewritten.

**`implanted_signal_profiler_scale` added.** Bracken: the read fraction.
MetaPhlAn implanted feature: `q_it / D_i` with `D_i = (1-F_i) + Q_i`.
Non-implanted: 0. A hard gate now asserts
`expected_abundance_profiler_scale = retained_baseline_profiler_scale +
implanted_signal_profiler_scale` for every row.
`implanted_signal_read_proportional` is the explicit alias of
`target_fraction_for_feature`, whose meaning is unchanged.

**Recovery numerator and denominator now share a scale.** Under
`profiler_scale`, recovery divides `response_signal_profiler_scale` by
`implanted_signal_profiler_scale`; it never compares a profiler-scale response
with a read fraction. Both are emitted as `response_signal_selected` and
`implanted_signal_selected`, with `reference_scale`, `selected_reference_type`
and the original row-level `reference_type`. A direct implanted target with a
missing, non-finite or non-positive implanted signal is a hard failure.

**R models no longer reject the profiler scale.**
`fit_continuous_dose_response.R` reads `implanted_signal_profiler_scale` and
`recovered_spike_signal_profiler_scale` from the endpoint table, validates the
profiler mapping, and uses the continuous implanted signal as the predictor.
Because `G_eff,i` varies by sample, the nominal categorical-dose nonlinearity
test is fitted only when a valid grouped dose variable is retained; otherwise it
is omitted and recorded in `reference_scale_manifest.tsv`. The read-proportional
model is behaviourally unchanged.
`fit_ideal_counterfactual_did.R` gains a required `--response-table` for the
profiler-scale path and takes `expected_abundance_profiler_scale` verbatim from
the migrated all-feature response table. It never recalculates `Q_i` in R.

**Target genome sizes.** `scripts/build_target_genome_sizes.py` measures the ten
lengths directly from the FASTAs referenced by `spikes/spike_panel.tsv`, with
per-FASTA SHA-256 and provenance, failing on a missing, empty, duplicated,
headerless or unexpectedly labelled entry. The FASTAs are not present locally, so
**the table has not been generated and is not audited**; André must run it.

**Tests.** 14 Python suites and the R DiD suite pass, including
`test_combined_profiler_reference.py` (12 cases covering the eight validation
rules, combined-profiler acceptance, the deliberately different recovery class,
Bracken and detection invariance and no-silent-fallback) and
`test_target_genome_sizes.py` (7 cases). Two vacuous or wrong assertions were
caught and fixed while writing them: the DiD residual merge originally matched
zero rows, and the analyzer's TSV path bypassed the reference remap.

**Not executed:** `fit_continuous_dose_response.R` — `mgcv` is absent, no image
is built and `apptainer` is unavailable. It was validated by `parse()` only, so
its profiler-scale path is unverified.

## 2026-09-18 — analyzer defect fixed; runner interfaces; checklist contract

The G_eff policy is unchanged. The fitted 3.10 Mb value remains unused anywhere.

**Defect: the analyzer rejected genuine profiler-scale builder output.**
Reproduced end to end as `[ERROR] line 2: inconsistent signed bounded error`.
Two causes, both now fixed:

1. The profiler-scale columns were remapped onto the canonical expected,
   retained and signal fields, but `signed_bounded_error` and
   `absolute_bounded_error` stayed read-proportional, so a genome-equivalent row
   carried an expected value from one scale and a bounded error from the other.
2. `analyze_perturbation_response.py` validated by reconstructing the retained
   baseline as `(1 - effective_total_fraction) * baseline_abundance_fraction`.
   That is the read-proportional formula and is wrong for MetaPhlAn, whose
   composition is renormalised by `D_i = (1 - F_i) + Q_i`.

**Correction.** `build_perturbation_response_input.py` now emits
`signed_bounded_error_profiler_scale` and
`absolute_bounded_error_profiler_scale`, validated finite and within [-1, 1] and
[0, 1] with `absolute == abs(signed)`. `reference_scale.py` maps both onto the
canonical names. The analyzer takes the supplied retained baseline as
authoritative and validates only scale-independent identities:
`expected = retained + implanted_signal`,
`response_signal = observed - retained`, `response_delta = observed - expected`,
`signed = (observed - expected)/(observed + expected)` or 0, and
`absolute = abs(signed)`. `effective_total_fraction` is retained as perturbation
metadata and no longer defines the selected retained abundance. The
read-proportional bounded-error columns are unchanged.

**Required regression.** `tests/test_builder_analyzer_profiler_scale.py` runs the
real builder and feeds its genuine Parquet to the real analyzer. It hand-checks
`q_it`, `Q_i`, `D_i`, retained baseline, implanted signal, expected abundance,
response signal and delta, and both bounded errors; asserts a non-implanted
MetaPhlAn feature uses the complete-mixture denominator; asserts Bracken and
detection are identical across scale selections; and asserts that a corrupted
profiler-scale bounded error makes analysis fail rather than pass silently.
A fixture that merely fabricates consistent analyzer input is not sufficient.

**Runner interfaces.** `run_yachida_definitive_analysis.sh` and
`run_crc_cohort_definitive_analysis.sh` now pass
`--reference-scale profiler_scale` to `fit_continuous_dose_response.R`; no
default was added inside the R script. `run_ideal_counterfactual_did.sh` requires
`REFERENCE_SCALE` (only `profiler_scale` or `read_proportional`), always passes
it, and under `profiler_scale` additionally requires a non-empty `RESPONSE_TABLE`
which it passes as `--response-table`. `read_proportional` does not require one.
The response table must come from the same cohort, population, assembly arm,
canonical input and endpoint generation as the DiD run. No shell DAG currently
invokes that runner, so no DAG stage needed adding; the documented callers were
updated. `tests/test_runner_reference_interfaces.py` covers all of this.

**Checklist contract.** C5 is `PENDING`, not a new status value. Its code is
implemented but definitive execution and sealing are outstanding, and that
distinction lives in the question, estimand, source, output and required-seal
columns. The allowed vocabulary was not expanded and
`test_manuscript_results_checklist.py` passes unchanged.

**Verification.** All 40 Python test files pass. `test_ideal_counterfactual_did.R`
passes. `test_continuous_model.R` **skips** — it prints
`[SKIP] mgcv unavailable`, so it did not execute and is not claimed as passing.
`fit_continuous_dose_response.R` parses. Every modified runner passes `bash -n`.

**Still cluster-dependent, and agents never run these:** real FASTA-derived
target genome sizes, real G_eff coverage, execution of both R model stages, and
definitive execution and sealing. Project decision documents must be updated
whenever a methodological or execution decision changes.

## 2026-09-20 — Frozen environment required for profiler-scale propagation

- **Decision:** use a new `ground_truth_analysis_v2_1.sif`; never mutate or
  overwrite the historical v1 or original-manuscript image.
- **Reason:** direct cluster checks established that both historical images lack
  Python DuckDB; the original manuscript image also lacks R `sandwich` and
  `arrow`. Bare login-node Python therefore cannot validate or execute the new
  response-table workflow.
- **Implementation:** versioned Apptainer definition and Conda environment;
  pinned Python 3.11, DuckDB 1.0.0, PyArrow 17.0.0, R 4.3.3, MaAsLin2 1.18.0,
  sandwich 3.1.1, arrow 17.0.0, mgcv and the prior manuscript dependencies. Build-time and
  post-install verification fail closed. The builder records the image SHA-256,
  source-file hashes, explicit Conda manifest, image inspection and R session
  information.
- **Status:** `OPEN` until André builds the image on the cluster and all four
  containerized regression tests pass. No profiler-scale propagation is
  authorized before then.

**Patch-level correction.** The first v2 image candidate passed its original
verifier but failed immediately when the real analyzer imported
`pyarrow.parquet`. R `arrow` and Python `pyarrow` are separate packages; only
the former had been declared. The incomplete image is retained as historical
evidence and is not patched or reused. Version 2.1 adds pinned PyArrow 17.0.0,
imports both `pyarrow` and `pyarrow.parquet` in the build verifier, and must be
built at a new path before rerunning the gates.

## 2026-09-20 — Remove historical-input discovery from G-eff propagation

- **Decision:** do not select `PRIOR_ENDPOINTS` or `ABUNDANCE_LONG` by searching
  historical `work/` trees.
- **Bracken invariant:** compare Bracken rows from the primary profiler-scale
  endpoints against the read-reference sensitivity endpoints generated in the
  same run from the same canonical input. Cohort, population, sample, target,
  dose and quantitative fields are compared exactly.
- **Abundance source:** deterministically rebuild the complete native abundance
  table from the checksum-locked Yachida canonical manifest, its exact
  `source_profile` files and the maintained alias table. Require its unit test,
  `SUCCESS`, exact schema and SHA-256 provenance before response construction.
- **Reason:** recursive NFS discovery is slow, historical endpoint tables may
  have incompatible cohort scope, and MapReduce abundance files are shards
  rather than guaranteed complete inputs. Same-run derivation removes all
  three ambiguities.
- **Status:** `DECIDED`; encoded in the cluster handoff.

**Interface correction after first controlled execution.** The response builder
requires the generated `biomarker_profile_manifest.tsv`, whose schema includes
`target_feature`; the raw canonical table intentionally does not contain that
derived field. The first development attempt correctly stopped at DuckDB bind
time after completing abundance expansion, both endpoint arms and the Bracken
check. The maintained runner now passes the generated manifest, asserts its
schema, and supports `RESUME_RUN_ROOT`. Completed abundance/endpoints are reused,
while an incomplete response-input directory is moved under `failed_attempts/`
before retry. No completed artifact is deleted or recomputed during this resume.
