# Next actions and execution order

**Purpose:** operational starting point for the next agent.
**Last updated:** 24 September 2026.

## Current operational checkpoint

- All three paper cohorts are complete and production-sealed. Yachida has 201
  samples, Feng 154, and Zeller 156. The Yachida assembly-sensitivity
  experiment is also sealed at 360/360 expected profiles.
- The topology-hardened CRC audits completed on 24 September 2026: Feng job
  `3097587` and Zeller job `3097588`, both `COMPLETED` with exit code `0:0` and
  empty error logs. Feng contains 154 baseline, 1,078 community, and 1,800
  independent profiles. Zeller contains 156 baseline, 1,092 community, and
  1,800 independent profiles. Every member of both production-seal checksum
  manifests verified `OK`.
- Upstream profiling must not be rerun. The next upstream-facing task is to
  build and review the checksummed evidence package specified in
  `UPSTREAM_EVIDENCE_PACKAGE_SPEC.md`; definitive downstream execution remains
  outstanding.
- The GUTBIOME Control/LR/HR profiling pilot is separate from the paper's
  three-cohort spike benchmark and must not enter definitive paper inputs.

Public-release documentation and deferred metadata/archive work are tracked in
`PUBLIC_RELEASE_CHECKLIST.md`. That checklist is the concise release-facing
task list; this file remains the detailed operational handoff.

**Provisional paper figure book (21 September 2026):** Six main figures and
seven supplementary candidates are inventoried in
`analysis_v2/PROVISIONAL_FIGURE_INVENTORY.tsv`; a checksum-tracked HTML review
book can be assembled with `scripts/assemble_provisional_figure_book.py`.
This is a layout and evidence-gap review, not a claim that all panels are ready.
The corrected three-cohort endpoint and recovery sources were run by André,
but their rendered previews remain provisional; Figures 3–6 still have
analysis/model gates, and the baseline-prevalence discrepancy is unresolved.

**Three-panel recoverability successor:** code and local fixtures now cover all
three cohorts and panels A–C; see `analysis_v2/THREE_COHORT_RECOVERABILITY_FIGURE.md`.
It is **not yet empirically run**. It requires a single corrected three-cohort
profiler-scale paired-endpoint table; the older mixed-provenance read-reference
table is not acceptable. Do not treat the existing two-cohort driver figure as
the updated result or bypass the new scale gate.

**Time-sensitive draft baseline figure:** `build_three_cohort_baseline_figure_input.py`
and `plot_three_cohort_baseline_discordance.R` now derive the four-taxon
prevalence/Wilson-interval and positive-only abundance panels directly from the
sealed three-cohort map-reduce `common_input` tables. They do not re-profile or
refit models. The Python/R fixture passed locally; the real three-cohort plot
has not yet been executed or reviewed and remains `DEVELOPMENT_ONLY`.

**Figure 6 in flight:** the sealed three-cohort development disease-model table
was located at `work/analysis_v2_three_cohort_mapreduce_dev_20260913_193010/models/disease/models/primary_disease_da_results.tsv` (932 MB; all three cohorts and both profilers). Before evaluating it, the community-dose axis was corrected: `spike_fraction_target` is the median per-member dose, while `spike_fraction_total` records the ten-member sum. Pull the corrected commit on the cluster before evaluation; use a compute allocation because the existing input is large. No real corrected Figure 6 result has been reviewed.

## Immediate result and next gate (21 September 2026)

The controlled Yachida `DEVELOPMENT_ONLY` propagation and residual-genome-size
audit completed successfully in
`work/geff_propagation_dev_20260920T003139Z`; every entry in
`run_checksums.sha256` verified. The next scientific gate is replication in
Feng and Zeller followed by the definitive three-cohort rebuild. Do not rerun
Yachida upstream profiling and do not treat the Yachida-only estimates as final
manuscript evidence.

### New: MetaPhlAn genome-size residual audit (20 September 2026)

`analysis_v2/scripts/audit_metaphlan_genome_size_residual.py` is implemented and
wired as resumable **stage 7** of `run_geff_propagation_development.sh`, after
the `reference_comparison/SUCCESS` gate. It answers one question: does the
apparent MetaPhlAn genome-size-dependent recovery bias disappear on the
genome-equivalent scale?

- Expected slope of the target-level median `log2(observed/expected)` against
  `log2(G_t)`: about **-1** before the correction (read-proportional
  sensitivity arm) and about **0** after it (genome-equivalent primary arm).
  Both expectations are written into the output and compared to the fit; the
  audit never assumes or forces either.
- The **implanted taxon is the statistical unit** — ten target-level points, not
  millions of observation rows. The bootstrap resamples **targets**, never
  observation rows, and is descriptive uncertainty over those ten taxa only
  (10,000 replicates, fixed seed 20260920).
- Genome sizes come only from the FASTA-measured `target_genome_sizes.tsv`. **No
  fitted 3.10 Mb constant is used anywhere**, and no genome size is estimated
  from a recovery outcome.
- **Corrected after Codex verification (20 September 2026).** Four fail-closed
  and provenance defects were fixed without changing the estimand, the
  target-level aggregation, the paired exclusion, the bootstrap, the seed, the
  scopes or the output names:
  1. the **selected** reference is now carried explicitly. Original row
     provenance (`reference_type`) and the selected estimand are different
     fields, so a MetaPhlAn sensitivity row that legitimately keeps
     `genome_equivalent` row provenance is now reported as the
     `read_proportional` estimand. The comparator emits
     `*_row_reference_type`, `*_selected_reference_type`, `*_reference_scale`
     and `*_selected_estimand`; the ambiguous `*_reference_type` names remain
     only as backward-compatible aliases and the audit does not read them;
  2. the audit **proves** its Yachida-only scope: every row's cohort is
     inspected and the distinct set must be exactly `{"yachida"}`. Blank,
     Feng-only, Zeller-only and mixed input all fail, and nothing is silently
     filtered;
  3. **every physical-identity component must be present and nonblank** before
     the duplicate-key check, so blanks can never collapse two distinct rows;
  4. **absolute relative error is validated** — numeric, finite, non-negative
     and equal to `abs(observed_over_expected - 1)` within
     `max(1e-12, 1e-9 * max(abs(actual), abs(expected), 1))` — and never
     recomputed and silently substituted. Paired-excluded rows are not rescued
     by it.
- Paired exclusion remains the accepted policy, unchanged.
- **Resume safety (20 September 2026).** The existing completed Yachida run's
  `reference_comparison` predates the explicit provenance columns. The runner no
  longer trusts a bare `SUCCESS`: `analysis_v2/lib/stage_compatibility.sh` +
  `scripts/check_stage_schema.py` verify the required files, exact
  tab-separated header columns and validation metrics before reuse, quarantine
  an incompatible directory under `$RUN_ROOT/failed_attempts/` and regenerate
  it. This applies to both `reference_comparison` and
  `metaphlan_genome_size_residual_audit`. The two target-recovery arms are
  inputs and are never regenerated or deleted, so **resuming that run will
  rebuild only the comparison and the audit.**
- **The real Yachida audit passed on 21 September 2026.** Across independent,
  community and pooled scopes, read-reference slopes were -0.9713, -0.9995 and
  -0.9962 (all bootstrap intervals excluded zero), whereas genome-equivalent
  slopes were 0.0132, 0.0068 and 0.0080 (all intervals included zero). All
  10,000 target-bootstrap replicates were valid in every fit. The audit used 10
  implanted taxa, proved the Yachida-only cohort scope, preserved Bracken
  identity, and used measured spike-FASTA lengths without a fitted genome-size
  constant or pseudocount.
- Status stays `DEVELOPMENT_ONLY` until Feng and Zeller replicate the finding
  and the three-cohort run validates it. The Yachida slopes may be quoted only
  with that explicit limitation.

Read `CLAUDE.md` or `CODEX.md` first; they are the same file. Record every
substantive decision in `analysis_v2/METHODS_DECISION_LOG.md` and update the
shared context in the same session.

## Immediate assignment

The profiler-scale implementation, Yachida propagation, paired reference
comparison, and residual-genome-size audit are complete and checksum-verified.
Next: verify the corrected Figure 6 path on real development data, execute the
real continuous model and DiD, and complete Feng/Zeller replication before
freezing publication figures or final Results.

**Progress, 17 September 2026.** The `G_j` source and `G_eff,i` derivation are
decided and validated against the production MetaPhlAn vJan25 database and all
201 sealed Yachida baselines; see §1. `scripts/derive_paired_endpoints.py`
implements the corrected target reference with additive fields and fail-closed
gates. The immediate blocker is now migration of the runners and downstream
consumers listed in §2, including the non-implanted reference.

Yachida is ready: 201/201 samples and 3,408/3,408 expected MetaPhlAn outputs
are verified. All ten exact spike FASTAs, lengths, and SHA-256 hashes are in
`CLAUDE.md`/`CODEX.md`. Zeller and Feng are not yet ready for the definitive
three-cohort run; always recheck their live state.

**Progress, 18 September 2026 — profiler-scale correction landed (local only).**
`analyze_perturbation_response.py` now honours the selected reference scale on
every quantitative path, not only in row validation:

- `ResponseRow.retained_baseline` returns the supplied selected retained
  baseline and never reconstructs `(1-F)o`;
- a deterministic pre-pass supplies one selected-scale implanted driver per
  independent observation (`q_it/D_i` for MetaPhlAn under `profiler_scale`),
  used by the operator design, the regression and the cross-cohort holdout
  prediction `retained_baseline + driver * slope`;
- `exact_superposition_rows(path, reference_scale)` selects one internally
  consistent column triple and fails closed on the wrong one.

Detection, Bracken and all read-proportional columns are unchanged, asserted by
`analysis_v2/tests/test_builder_analyzer_profiler_scale.py` (17 tests, hand-
computed numbers). The stale CLAUDE.md line saying the downstream consumers were
"not yet converted" has been corrected, as has the obsolete header comment in
`fit_continuous_dose_response.R` claiming `profiler_scale` fails closed.

**Still open before any definitive run:**

1. **DONE (19 September 2026):** Codex re-verified the profiler-scale
   correction and the exactly-one-driver fail-closed guard.
2. `mgcv` is absent locally, so `fit_continuous_dose_response.R` has been
   parse-checked only and has **never been executed**. Do not report otherwise.
3. The ten-target genome-size table is still ungenerated (cluster-only FASTAs).
4. **IMPLEMENTED, awaiting real execution — Figure 6 community correction.**
   Community rows collapse to one physical mixture per dose, all panel taxa are
   excluded after alias canonicalization, and missing members or inconsistent
   repeated fits fail closed. Existing Figure 6 output remains unlicensed until
   the corrected evaluator and plot are rerun and reviewed.

## 0. Protect and understand active work

Read, in order:

1. `CLAUDE.md` or `CODEX.md`;
2. this file;
3. `WORK_HANDOFF.md`;
4. `PREPRINT_REVISION_ROADMAP.md`;
5. `analysis_v2/STATISTICAL_ANALYSIS_PLAN.md`;
6. `analysis_v2/METHODS_DECISION_LOG.md`;
7. `analysis_v2/PROFILER_SEMANTICS.md`;
8. `analysis_v2/ENDPOINTS.md`.

Run:

```bash
git status --short --branch
git diff --stat
git ls-files --others --exclude-standard | sort
```

Preserve all current modified and untracked work. Do not reset, clean,
overwrite, or make a broad commit without auditing provenance and tests.
Confirm `CODEX.md` still resolves to `CLAUDE.md`.

## 1. Effective-community-genome-size policy — completed

**Status: DECIDED and validated on sealed Yachida.** `G_eff` is sample-wide
(one row per `(cohort, sample_id)`, blank `analysis_population` for the generic
downstream lookup), uses exact full-lineage terminal-SGB matching, and is
derived only from unspiked MetaPhlAn baselines. The production vJan25 database
contained 58,216 usable terminal-SGB genome lengths. All 201 baselines had SGB
rows; exact-lineage mapping coverage was 1.0 for every sample, with zero
unmapped features and zero exclusions. Both the prespecified 95% primary and
90% sensitivity gates passed. `G_eff` ranged 2.630--4.844 Mb (median 3.477 Mb).

Audit bundle:
`work/yachida_geff_audit_20260917T144215Z`. Database SHA-256:
`e7d23a73a7959b4f41af0bbe403f4b5bbb7c1879528d376d146ee5294515df9a`.
André ran the cluster audit manually; agents did not access the cluster.

The estimator, fail-closed policy, diagnostics, provenance, reproduction
commands, and limitations are in
`analysis_v2/METAPHLAN_EFFECTIVE_GENOME_SIZE.md`. This completes derivation of
the reference input only; it does **not** by itself license final MetaPhlAn
recovery claims.

The frozen policy specifies:

1. authoritative genome-size source for baseline MetaPhlAn taxa/SGBs;
2. feature-to-genome-size identifier mapping;
3. sample-specific, cohort-specific, or global scope;
4. mathematical weighting rule;
5. minimum mapped-abundance coverage;
6. treatment of unmapped taxa and non-species/unclassified mass;
7. diagnostics and failure thresholds;
8. prespecified sensitivity alternatives;
9. proof that post-spike outcomes are not used to choose `G_eff,i`.

Frozen sample-specific estimator:

`G_eff,i = sum_j(c_ij * G_j) / sum_j(c_ij)`

where `c_ij` is baseline MetaPhlAn genome-equivalent-like relative abundance
and `G_j` is the production MetaPhlAn database representative genome size.

Do not fit `G_eff` to make recovery unbiased, use the provisional 3.10 Mb fitted
value as truth, substitute `f/mapped_fraction`, or call whole-genome scaling an
exact reproduction of marker normalization.

Completed deliverables:

- a dated `DECIDED` entry in `METHODS_DECISION_LOG.md`;
- a module document with formula, inputs, exclusions, diagnostics, and
  sensitivities;
- a versioned genome-size mapping table with provenance;
- a mapping-coverage audit on sealed Yachida baselines.

## 1b. Propagate the decided reference scale (partially done)

**Done:** four runners migrated through `analysis_v2/lib/metaphlan_reference.sh`
(fail-closed, no default, separate primary/sensitivity directories, provenance).
Cluster commands ready in `analysis_v2/CLUSTER_HANDOFF_GEFF_PROPAGATION.md` —
André runs them.

**Done:** non-implanted MetaPhlAn reference, `summarize_target_recovery.py`,
`analyze_perturbation_response.py`, runner hardening, repaired baseline test,
corrected cluster handoff. **Remaining:** `fit_ideal_counterfactual_did.R`
profiler_scale needs a per-member implanted-fraction input; report/plot scripts
still inherit whatever scale their input was built with.

**18 Sep 2026.** Analyzer profiler-scale defect fixed and regression-tested
(real builder -> real analyzer). Runner interfaces migrated; C5 set to PENDING.
All 40 Python test files pass; `test_ideal_counterfactual_did.R` passes;
`test_continuous_model.R` SKIPS (no mgcv).

**Next action:** Codex re-verification. Then André, personally: run
`analysis_v2/scripts/build_target_genome_sizes.py` on the cluster, run the
continuous-model contract test inside the analysis image, then set
`PRIOR_ENDPOINTS`, `ABUNDANCE_LONG` and `TARGET_GENOME_SIZES` and follow
`analysis_v2/CLUSTER_HANDOFF_GEFF_PROPAGATION.md`. Agents never run cluster
commands.

**Untested:** `fit_continuous_dose_response.R` profiler-scale path (`mgcv`
absent locally).

## 2. Implement profiler-specific endpoints

**Status:** `scripts/derive_paired_endpoints.py` is **done** (additive fields,
`reference_type`, `--metaphlan-reference` required with no default, Q_i grouped
per `profile_id`, community-completeness gate, tests in
`tests/test_paired_endpoints_profiler_scale.py`). The remaining scripts below
are **not started**, and neither is the non-implanted reference `e_MP,ij`.

Read `analysis_v2/INPUT_CONTRACT.md`, `analysis_v2/ENDPOINTS.md`,
`analysis_v2/CONTINUOUS_MODEL.md`,
`analysis_v2/PERTURBATION_RESPONSE_ATLAS.md`, and
`analysis_v2/IDEAL_COUNTERFACTUAL_DID.md` before editing.

Preserve old endpoints under explicit `*_read_reference` names. Add fields for
reference type, target genome size, effective community genome size, implanted
genome-equivalent fractions, profiler-scale expected abundance, retained
baseline, and residual.

Bracken primary reference:

`e = (1-F)o + f`.

MetaPhlAn first-order reference:

`q_it = f_it * G_eff,i / G_t`

`Q_i = sum_k(f_ik * G_eff,i / G_k)`

Target:

`e_MP,it = ((1-F_i)o_it + q_it) / ((1-F_i) + Q_i)`

Non-implanted feature:

`e_MP,ij = (1-F_i)o_ij / ((1-F_i) + Q_i)`.

Likely affected code:

- `analysis_v2/scripts/derive_paired_endpoints.py`;
- `analysis_v2/scripts/build_perturbation_response_input.py`;
- `analysis_v2/scripts/analyze_perturbation_response.py`;
- `analysis_v2/scripts/summarize_target_recovery.py`;
- `analysis_v2/scripts/fit_continuous_dose_response.R`;
- `analysis_v2/scripts/fit_ideal_counterfactual_did.R`;
- dependent reports and plots.

Add tests for equal and unequal genome sizes, individual/community mixtures,
zero baseline, non-implanted features, incomplete mappings, unchanged Bracken,
distinct reference outputs, and fail-closed genome metadata.

## 3. Validate on sealed Yachida

Run in `DEVELOPMENT_ONLY` mode using a new output directory.

Acceptance criteria:

1. 201 samples and all 3,408 MetaPhlAn outputs are represented.
2. Target FASTA hashes match the audited values.
3. Every included sample has `G_eff,i` or a controlled exclusion.
4. Bracken endpoints and all detection results are unchanged.
5. Read-reference and profiler-scale results cannot be confused.
6. Genome-size dependence is retested at sample/cohort level with uncertainty.
7. Prespecified `G_eff` sensitivities are reported.
8. Diagnostics, source tables, provenance, checksums, session information, and
   `SUCCESS` are produced.

If residual genome-size dependence remains, investigate marker representation,
mappability, strain divergence, copy number, and database composition. Do not
tune `G_eff` until it disappears.

## 4. Finish and seal Feng and Zeller

This can proceed operationally in parallel, but it must not influence selection
of the MetaPhlAn correction.

Required final counts:

- Yachida: 201/201, already sealed;
- Zeller: 156/156;
- Feng: 154/154.

Read `datasets/CRC_ROLLING_PRODUCTION.md`,
`datasets/CRC_ROLLING_SUPERVISOR.md`, `datasets/README.md`,
`datasets/yachida/README.md`, and the definitive runbooks under `analysis_v2/`.

Require `.verified`, retained-output receipt, checksums, cohort audit, dataset
seal, and nonempty `SUCCESS`. Slurm `COMPLETED` is insufficient.

## 5. Freeze identity and profiler semantics

Freeze `analysis_v2/feature_equivalence_aliases.tsv`, preserve raw profiler
identifiers, prohibit broad prefix collapsing, keep
`Fusobacterium nucleatum_E` separate unless independently justified, and audit
native fields/database versions. Read `analysis_v2/PROFILER_SEMANTICS.md` and
the corresponding decision-log entries.

## 6. Build the definitive three-cohort input

Use sealed strict-production inputs only. Validate exact cohort/condition
counts, `Control`–`Adenoma`–`CRC` order, population separation, pairing,
fractions, complete panels, raw/canonical identities, native units, paths,
receipts, checksums, and zero unexplained missingness.

Read `analysis_v2/INPUT_CONTRACT.md`,
`analysis_v2/THREE_COHORT_DEVELOPMENT.md`, and `REPRODUCING.md`.

## 7. Run the frozen definitive analysis once

In `DEFINITIVE` mode run:

1. detection;
2. profiler-specific quantitative response;
3. artificial biomarker recovery;
4. corrected non-implanted response atlas;
5. corrected community superposition;
6. native CRC biomarker robustness;
7. corrected ideal-counterfactual distortion;
8. cross-cohort directional replication;
9. required sensitivity analyses;
10. cross-cohort synthesis and publication report.

Read `analysis_v2/STATISTICAL_ANALYSIS_PLAN.md`,
`analysis_v2/ANALYSIS_POLICY.tsv`, module documents, and
`analysis_v2/REPORTING.md` first.

## 8. Generalization and falsification analyses

If feasible and frozen prospectively:

- leave-one-cohort-out response transfer;
- cross-cohort community reconstruction;
- phenotype-permutation truth-anchored fragility;
- perturbation stability versus external replication;
- reference/database representation sensitivity;
- sequencing-depth sensitivity.

Retain negative results.

## 9. Final figures and manuscript

Proposed main figures:

1. design and estimands;
2. detection and corrected quantitative response;
3. taxonomic cross-talk;
4. independent-to-community reconstruction;
5. CRC biomarker stability;
6. cross-cohort transfer and perturbation stability.

Use `PREPRINT_REVISION_ROADMAP.md` for the paper structure and claim-to-evidence
plan. Generate every final product from the same sealed definitive bundle.

## 10. Submission audits

Perform separate statistical, taxonomic-identity, reproducibility/provenance,
and claim-to-evidence audits. Every final figure needs source data, exact
denominators, effect sizes and uncertainty, status, provenance, checksums, and
a regeneration command.

## Completion rule

The project is submission-ready only after the MetaPhlAn correction is
validated, all cohorts are sealed, the definitive workflow and sensitivities
pass, the manuscript is regenerated from definitive evidence, and all four
submission audits pass.
