# Next actions and execution order

**Purpose:** operational starting point for the next agent.
**Last updated:** 19 September 2026.

## Immediate infrastructure gate (20 September 2026)

The profiler-scale code is locally verified, but neither historical cluster
image contains the complete new dependency set. Build and verify
`ground_truth_analysis_v2.sif` using
`build_ground_truth_analysis_v2_container.sh`, then run all four containerized
regressions in `analysis_v2/CLUSTER_HANDOFF_GEFF_PROPAGATION.md`. Do not resume
the `DEVELOPMENT_ONLY` propagation until every image and regression gate passes.
Do not overwrite the historical v1 or original-manuscript images.

Read `CLAUDE.md` or `CODEX.md` first; they are the same file. Record every
substantive decision in `analysis_v2/METHODS_DECISION_LOG.md` and update the
shared context in the same session.

## Immediate assignment

The profiler-scale implementation and maintained runners are migrated and
Codex-verified. The next task is to prepare a scoped commit and push, then have
André personally perform the controlled `DEVELOPMENT_ONLY` execution in
`analysis_v2/CLUSTER_HANDOFF_GEFF_PROPAGATION.md`. Do not rerun upstream
profiling or regenerate poster figures.

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
4. **PENDING, not started — Figure 6 implanted-taxon exclusion.** Community
   profiles exclude only the focal target instead of all ten implanted taxa;
   community profiles may be counted repeatedly across `target_label` rows; and
   feature aliases must be canonicalized before target exclusion. Recorded in
   `analysis_v2/METHODS_DECISION_LOG.md` (2026-09-18) as `OPEN`. This was
   deliberately excluded from the reference-scale task and needs its own task.

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
