# Downstream scientific scope and readiness audit

Date: 29 September 2026. Reviewed code: `dbf3d6a`.
Status: review and proposed execution priorities, not an amendment to the frozen
analysis policy. No statistical methods or results were changed for this audit.
Do not submit a full downstream run until the execution blockers below are fixed
and the intended paper modules are agreed.

## Evidence reviewed

- Local Overleaf `sn-article.tex`, titled *Ground-truth in silico spike-ins reveal
  limits of microbiome biomarker recovery in colorectal cancer*, especially its
  abstract, Results, and Methods. `new_version.tex` is an older/different draft
  and was not treated as the authoritative manuscript. The manuscript stays
  outside the code repository.
- `PREPRINT_REVISION_ROADMAP.md`, `NEXT_ACTIONS.md`, project context,
  `STATISTICAL_ANALYSIS_PLAN.md`, `ANALYSIS_POLICY.tsv`, model specifications,
  `METHODS_DECISION_LOG.md`, `MANUSCRIPT_RESULTS_CHECKLIST.tsv`, and the provisional
  figure inventory.
- The shared cohort runner, MetaPhlAn reference helper, quantitative model,
  calibration-linkage report, atlas, ideal-counterfactual, calibration,
  reliability, synthesis, and publication entry points.
- User-provided cluster results: unified upstream seals, evidence package, and
  shared preflights 3097753–3097755. No agent accessed the cluster.

## 1. What is actually ready

Upstream production is complete: 201 Yachida, 154 Feng, and 156 Zeller samples
(511 total). The additional assembly-sensitivity experiment has 360 profiles.
The three upstream seals use a common interface and the upstream evidence
package passed checksum verification. Bracken's recorded 100-base setting and
threshold 10 were verified. Extra parameter files were quarantined attempts,
not missing production profiles.

All three shared preflights completed with exit 0 and empty error logs:

| Cohort | Job | Biological samples | Canonical rows |
|---|---:|---:|---:|
| Yachida | 3097753 | 201 | 36,360 |
| Feng | 3097754 | 154 | 28,840 |
| Zeller | 3097755 | 156 | 29,160 |

These validate the common manifest handoff and 94,360 canonical rows. They do
not validate the model stages, corrected reference-table availability for every
cohort, or completion of the revised paper's entire analysis plan.

## 2. What needs to change from the preprint

### Quantitative recovery must use the declared profiler reference

The old read-based expectation is not the primary quantitative reference for
MetaPhlAn in the revised policy. The genome-equivalent correction affects
targets, retained background, non-implanted responses, community reconstruction,
and ideal-counterfactual distortion. It is a model of the expected measurement,
not direct cellular ground truth.

Documented Yachida development evidence illustrates the consequence: community
good recovery changed from 5.38% to 34.01%, and independent good recovery from
5.83% to 34.61%. These figures remain Yachida-only development evidence. They
cannot become final numerical claims without the three-cohort rebuild. The
preprint's strong MetaPhlAn quantitative-recovery ranking must be reassessed.
Detection itself is a separate endpoint and is not changed by that correction.

### The paired design must be respected

The preprint explicitly used unpaired MaAsLin2 spike-status models even though
the same biological samples supplied both states. The revised artificial-marker
model uses within-sample log2 differences, a fixed pseudocount, a frozen feature
family, and feature-wise tests over biological samples. Disease association is
a different model with age/sex adjustment and HC3 uncertainty.

Switching software back to MaAsLin2 is not a prerequisite. The paired scientific
contrast, experimental unit, feature universe, assumptions, diagnostics, and
multiplicity policy are what must be justified and consistently implemented.
The existing paired model should be checked before any replacement is proposed.

### Filtering claims need stricter interpretation

The preprint's artefact exclusion removed entries from previously fitted
significant-call lists, without refitting MaAsLin2 or BH. It protected known
implanted targets, so unchanged target recall was guaranteed by construction.
Its reciprocal cross-study application provides some transfer evidence, but
does not establish a general feature-selection procedure that preserves unknown
biomarkers. The revision's blind development screen was largely negative or
asymmetric. A useful outcome is to document those limits.

The preprint's calibratability analysis already included grouped sample
cross-validation and reciprocal cohort transfer. Preserve that useful design
principle, but regenerate quantitative conclusions against the corrected
reference. Do not characterize every historical analysis as unvalidated.

### Separate three meanings of "biomarker"

1. A taxon is reported by a profiler: operational detection.
2. A known spike is significantly enriched relative to its paired baseline:
   artificial-marker recovery.
3. A taxon is associated with CRC or adenoma relative to Control: native
   disease association.

Native biomarker losses/gains under perturbation describe stability, not known
biological truth. A changed q-value can arise through effect, uncertainty, and
the multiplicity family. Include continuous effects and an ideal-reference
comparison, rather than using call transitions alone as evidence of distortion.

## 3. Recommended scientific questions

### A. Can the workflows detect and quantify controlled additions?

Primary: unconditional paired quantitative response on the profiler-specific
scale, and native detection across dose. Preserve zero reports and negative
baseline-adjusted responses. Report cohort-specific effects and intervals.
Show taxon-specific patterns and clinical background without treating repeated
dose rows as independent samples.

Native nonzero abundance can reflect a taxon already present at baseline.
Therefore baseline detection should be visible, and a baseline-negative
descriptive stratum would help distinguish first detection from incremental
recovery. Adding such a stratum is a proposal, not an already frozen analysis.

### B. Is the implanted signal recovered as an artificial biomarker?

Use the paired spiked-versus-baseline model, with pooled clinical backgrounds
primary and stratified results secondary. Report effects, uncertainty, target
recall, enriched non-target burden and denominators. Use the frozen primary
q <= 0.05 and sensitivity q <= 0.10. Minimum tested significant dose is an
empirical statistic, not a universal analytical detection limit. It depends on
sample count, baseline, feature family, model, and dose grid.

### C. What happens to non-implanted features, and do mixtures behave additively?

Measure responses relative to the profiler-specific retained baseline. Compare
observed community profiles with composition-only expectation and reconstruction
from independent additions. Use matched samples and supported doses for a fair
independent/community comparison. Same-sample reconstruction must be labelled
as such; held-out cohort prediction is a stronger and separate question.

### D. Do these responses change actual disease-marker inference?

CRC versus Control is primary; Adenoma versus Control is secondary. Keep
full-community and independent-subset baselines separate. Report native effect
changes, retention, loss, gain, and direction changes. Exclude all directly
implanted community members from non-target summaries and count each physical
mixture once. Fit the ideal-counterfactual distortion model to distinguish
departures from expected compositional response from q-threshold crossings.

### E. Which findings generalize across cohorts?

Show cohort-specific estimates, directional external replication and, where
specified, meta-analysis. These answer different questions. Any link between
perturbation stability and external replication is an analysis to test, not an
assumed benefit. With three cohorts, report the individual transfer directions
and uncertainty rather than promising universal transportability.

Filtering, abundance correction, an automatic reliability cutoff, and a
null-label test of whether fragility identifies false discoveries are optional
extensions. The core paper need not succeed at correction to be valuable.
They become necessary if the paper claims improved inference or false-positive
identification through a learned correction/filter.

## 4. Development versus execution inventory

| Component | Implementation/evidence found | Remaining work |
|---|---|---|
| Upstream and shared canonical input | Real three-cohort seals and preflights passed | Preserve evidence; no upstream rerun |
| Baseline discordance and taxon identity | Identity lock and draft-figure code | Rebuild on final exact IDs; reconcile historical prevalence plots |
| Detection models | Specified, code and fixture tests; in common runner | Real fits, convergence/separation checks, scientific interpretation |
| Profiler-specific quantitative recovery | Corrected code and Yachida development validation | All-cohort reference tables, real fits, scale comparison and diagnostics |
| Paired artificial biomarkers | Specified and implemented; common runner | Real final fits, feature-family/denominator and q-value review |
| Non-implanted response atlas | Development code/results | Correct launcher arguments, integrate selected reference and final execution |
| Independent/community reconstruction | Implemented within response analysis | Integrate final matched inputs, uncertainty, explicit reconstruction/transfer labels |
| Native disease-marker propagation | Implemented; common runner | Final covariate/sample panels, zero unexplained exclusions, corrected community summaries |
| Ideal-reference distortion | Model and runner exist | Integrate with same-scale response table and execute on final data |
| Cross-cohort disease meta-analysis | Separate synthesis entry point exists | Final cohort results, checksum verification and diagnostics |
| Directional replication and reliability association | Development analyses exist | Final-data rebuild; reliability runner explicitly blocks definitive mode |
| Filtering and abundance correction | Experimental code and negative/asymmetric screens | Decide scope; calibration runner explicitly allows development only; held-out refits needed for improvement claims |
| Sensitivities and publication graphics | Multiple scripts and provisional assets | Pick required set, wire it, execute, review, unify style and seal source data |

"Implemented" and a synthetic-test PASS do not mean "run and validated on all
three final cohorts." No full final-cohort modelling was reported in the latest
cluster outputs.

## 5. Concrete blockers found in the present execution path

1. **Endpoint directory conflict.** `run_cohort_definitive_analysis.sh:56`
   creates `endpoints`; at line 89 it calls `derive_endpoints_with_references`.
   The helper in `lib/metaphlan_reference.sh:112` rejects an existing endpoint
   directory. The preflight exits before reaching it. The assembly-sensitivity
   runner has the same directory-creation pattern. Repair and exercise the
   actual caller/helper interface before full jobs.
2. **Reference prerequisites are outside preflight.** The helper requires
   `TARGET_GENOME_SIZES` and `EFFECTIVE_GENOME_SIZES`, including sample coverage.
   The three successful preflights did not reach that gate. Yachida reference
   derivation is documented; complete final Feng/Zeller coverage is not proved
   by these preflights. Audit/provide those tables and move the relevant checks
   ahead of model fitting.
3. **One report still selects the old reference.**
   `scripts/link_calibration_to_biomarkers.R` reads `response_ratio`,
   `recovered_spike_signal`, and the old reference-error fields. The corrected
   endpoint schema preserves those fields as read-reference quantities. This
   report is called by the common runner without a scale choice. It must use
   explicit selected-scale fields or be deliberately labelled sensitivity-only.
4. **The main manuscript plan exceeds the common runner.** The common runner
   does not invoke the response atlas, community reconstruction,
   ideal-counterfactual model, reliability/replication analysis, or genome-size
   residual audit. The subsequent synthesis script primarily combines disease
   results. A successful cohort package would therefore still leave major
   proposed manuscript questions unanswered.
5. **The generic atlas launcher is stale.**
   `run_perturbation_response_atlas.sh` does not forward the reference selection
   and genome-table arguments now needed by the builder/analyzer/summary chain.
   The corrected development propagation runner is evidence that the component
   chain can run, but is not a substitute for integrating a final workflow.
6. **Scientific documentation and evidence labels have drifted.** Several
   model documents still define only `(1-F)o + f`, despite corrected code;
   task notes contain superseded upstream and execution states. The claim
   checklist uses IMPLEMENTED for modules that are not yet final evidence.
   Synchronize the authoritative plan and stage inventory before publication.

These findings come from source inspection. This audit has not fixed these
issues, run new statistical models, or inspected unprovided cluster outputs.

## 6. Decisions to review before final execution

- Keep the existing paired/artificial and adjusted/native contrasts distinct.
  Specify the experimental unit, feature family, effect, interval and q-value
  family for each module in one short analysis matrix.
- The current quantitative/detection GAMs have pooled profiler effects, a
  condition main effect, and sample/target random intercepts. These alone do
  not estimate every taxon- or background-specific interaction implied by a
  broad narrative. Decide which heterogeneity claims are descriptive and
  whether additional prespecified contrasts are required. Review residual and
  within-sample dependence diagnostics before treating precision as settled.
- Preserve nominal dose/rank for grouping and achieved/profiler-scale dose for
  quantitative calculations. Under corrected scale the continuous model can
  omit its categorical nonlinearity check if nominal grouping is unavailable;
  resolve that input contract explicitly.
- Keep primary/sensitivity choices fixed. Document any change after earlier
  development-data inspection honestly as a dated revision. The revised study
  has substantial development evidence already; do not describe it as untouched
  prospective validation.
- Distinguish per-context BH control from a global claim across every context.
  A global distortion/fragility claim requires an appropriate aggregate or
  hierarchical analysis; do not pool all significant contexts as independent
  discoveries.
- Require held-out training and full model/BH refitting if filtering/correction
  is retained as an improvement claim. Do not optimize a rule on the cohort
  used to evaluate it.

## 7. Recommended work order and figures

1. Agree the five scientific questions above and mark correction/reliability
   extensions optional. Reconcile the frozen policy, model documents and code.
2. Fix the six integration/documentation gaps, add meaningful stage-interface
   tests, and extend preflight to verify all reference inputs.
3. Build/verify reference coverage and corrected endpoints for all three
   cohorts. Reproduce the Yachida scale comparison in Feng and Zeller.
4. Run one representative real model stage with full diagnostics. Resolve
   prespecified failures transparently, then execute the agreed core pipeline
   for all cohorts with a common stage inventory.
5. Run cross-cohort synthesis, mandatory sensitivities, and publication
   reporting. Review final source tables before writing numerical conclusions.
6. Seal the complete downstream evidence package for the paper and Zenodo.

Proposed main-figure architecture, subject to the scope decision:

1. design, sample flow, and distinct measurement references;
2. detection and corrected quantitative recovery;
3. paired artificial-marker recovery and non-target discovery burden;
4. non-implanted responses and independent-to-community reconstruction;
5. native disease-marker stability with ideal-reference distortion;
6. cross-cohort consistency/replication, if it adds a clear final result.

Baseline case studies, recovery classes, minimum significant dose, all-taxon
panels, reference-scale/assembly sensitivities, diagnostics, and exploratory
correction belong mainly in the supplement. Use one versioned plotting theme,
fixed profiler/condition colours and taxon order, comparable axes, stated
denominators and intervals, and a source-data file for every exported figure.
Do not substitute visual consistency for scientific validation.

The upstream Methods can be drafted now. Downstream Methods should follow the
agreed analysis matrix; final numerical Results and abstract claims should wait
for the complete validated outputs.
