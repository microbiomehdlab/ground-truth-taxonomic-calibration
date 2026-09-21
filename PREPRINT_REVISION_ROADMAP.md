# Revision roadmap for a rigorous journal submission

**Working manuscript:** *Ground-truth in silico spike-ins reveal limits of microbiome biomarker recovery in colorectal cancer*
**Target level:** Nature Communications or a journal with comparable expectations
**Document status:** strategic revision plan, updated 21 September 2026
**Source manuscript reviewed:** bioRxiv version 1, DOI 10.64898/2026.08.02.742082
**Important evidence warning:** unless explicitly stated otherwise, current three-cohort downstream results are `DEVELOPMENT_ONLY`. They combine definitive Yachida output with historical Feng and Zeller output. Final estimates and figures must be regenerated from sealed strict-production outputs for all three cohorts.

### Evidence update: profiler-scale correction

The sealed Yachida development propagation established that the abundance-scale
choice materially affects MetaPhlAn recovery. Relative to the preserved
read-proportional sensitivity, the genome-equivalent primary analysis increased
good recovery from 5.38% to 34.01% for community perturbations and from 5.83%
to 34.61% for independent perturbations. Median absolute relative error fell
from 66.49% to 13.43% and from 57.67% to 12.32%, respectively. Bracken was
identical in all 15,870 paired Bracken observations.

Accordingly, the preprint-level claim that essentially no MetaPhlAn targets
recover well must not be retained as a primary conclusion. The revised paper
should state, provisionally, that aligning ground truth with the profiler's
genome-equivalent-like abundance scale removes a large systematic component of
apparent error, but meaningful taxon- and context-dependent residual error
remains. These are Yachida-only `DEVELOPMENT_ONLY` results and require
Feng/Zeller replication before final prose or figures are frozen.

The checksum-verified Yachida residual-genome-size audit independently supports
that interpretation. The slope of target-level median
`log2(observed/expected)` against `log2(target genome size)` was approximately
-1 under the read-proportional sensitivity (-0.971 independent, -0.999
community, -0.996 pooled), but approximately 0 under the genome-equivalent
primary reference (0.013, 0.0068, 0.0080). Read-reference bootstrap intervals
excluded zero and all corrected intervals included zero. This is compelling
Yachida-only validation of the scale correction, not yet a three-cohort claim.

Paper drafting may proceed for the Introduction, experimental design, frozen
methods, and manuscript structure. Final Results, abstract claims, Discussion,
and publication figures remain blocked on cross-cohort replication, stratified
sensitivity review, real model execution, corrected Figure 6 development
execution, and the definitive three-cohort rebuild.

## 1. Executive recommendation

The revised paper should not be presented as the old two-cohort preprint with one additional cohort. The scientific advance is now broader:

> Controlled read perturbations can expose profiler-specific detection limits, quantitative bias, taxonomic cross-talk, mixture non-additivity, and instability of downstream CRC biomarker conclusions that cannot be diagnosed from ordinary case-control benchmarking alone.

This framing is stronger because it connects the analytical ground truth to a general measurement problem. The spike-ins are not merely artificial biomarkers to be recovered. They are controlled interventions that let us map how a taxonomic profiling system responds and how those responses propagate into downstream inference.

The paper should make three claims, in this order:

1. **Measurement response is profiler-, taxon-, dose-, and background-dependent.** Detection alone is insufficient; quantitative response must be measured continuously.
2. **Controlled perturbations reveal structured collateral responses and mixture behavior.** Non-implanted taxa can move systematically, and independently measured response components do not always reconstruct community perturbations equally well.
3. **These measurement properties affect CRC biomarker inference, but they do not themselves establish biological truth.** Biomarkers can be stable or perturbation-sensitive, and transfer between cohorts is limited. Stability and external replication are complementary evidence dimensions.

The manuscript must avoid three overclaims:

- added read fraction is not cellular abundance, biomass, extraction efficiency, or colonization;
- a changing non-implanted taxon is not automatically a biological false positive;
- a biomarker that disappears after perturbation is perturbation-sensitive, not proven false.

## 2. What has changed since the preprint

### 2.1 Scope and data

| Component | Preprint | Revised project |
|---|---|---|
| Cohorts | Feng and Zeller | Feng, Yachida, and Zeller |
| Clinical metagenomes | 310 | 511 when all definitive cohorts are complete |
| Primary perturbation types | individual and community | retained, but explicitly analyzed as distinct experimental populations |
| Inference | dominated by recovery summaries and several unpaired comparisons | paired, baseline-adjusted analysis with biological sample as the inferential unit |
| Feature scope | implanted taxa and selected apparent artefacts | implanted taxa plus the complete eligible profiler-reported feature universe |
| Downstream question | recovery, filtering, and calibration | measurement response, cross-talk, mixture reconstruction, biomarker robustness, and transfer |

### 2.2 Conceptual improvements already made

The revised work has already corrected or expanded several important aspects of the preprint:

1. Added the Yachida cohort, providing a third biological background and a better basis for assessing generalization.
2. Defined the truth correctly as implanted paired-read identity, pair count, and final-library read fraction.
3. Replaced simple spiked-versus-group comparisons with matched baseline-adjusted endpoints.
4. Treated the biological sample, rather than profile rows or dose observations, as the replication unit.
5. Separated the full community population from the frozen independent-spike subset.
6. Made continuous quantitative response primary and threshold recovery classes secondary.
7. Added native detection models rather than treating zero and nonzero values only descriptively.
8. Added a target-to-non-implanted response atlas to measure taxonomic cross-talk.
9. Added an additive community-reconstruction analysis comparing single-target response information with a composition-only reference.
10. Added native CRC biomarker retention, loss, gain, effect fidelity, direction stability, and induced-call analyses.
11. Added cross-cohort directional replication as an independent evidence annotation.
12. Added an ideal read-proportional counterfactual and a phenotype-dependent distortion analysis.
13. Added explicit development/definitive statuses, input contracts, manifests, diagnostics, checksums, and fail-closed completion gates.
14. Corrected taxon-identity handling after discovering that broad alias collapsing could inflate *F. nucleatum* prevalence.
15. Rejected the simple conclusion that a reusable static taxon blacklist is broadly supported; the development-only cross-cohort filtering screen was mostly negative or asymmetric.

These changes are scientifically meaningful. They should be described as a redesigned analytical framework, not as cosmetic updates.

## 3. What is complete, provisional, and missing

### 3.1 Implemented and methodologically specified

The following components have code, written analysis policies, or development outputs:

- canonical paired endpoint derivation;
- expected read-proportional abundance and baseline-retained abundance;
- continuous recovery and native detection models;
- independent and community perturbation separation;
- paired artificial-biomarker models;
- non-implanted response atlas;
- same-sample community superposition analysis;
- native disease-biomarker propagation and transition ledgers;
- cross-cohort biomarker replication matrices;
- perturbation-reliability summaries;
- ideal-counterfactual difference-in-differences model;
- experimental cross-cohort abundance calibration;
- assembly-choice sensitivity machinery;
- profiler-semantics audits and denominator-harmonized sensitivity analyses;
- strict upstream validation, receipts, checksums, and completion seals;
- a prespecified statistical analysis plan and methods decision log.

### 3.2 Upstream evidence status

- **Yachida:** strict production complete and sealed. It includes 201 baselines, 1,407 community profiles, and 1,800 independent profiles.
- **Zeller:** strict production was nearly complete at the latest recorded checkpoint; current cluster state must be audited and the cohort sealed.
- **Feng:** strict production has been submitted through the rolling scheduler. It is not yet a sealed definitive input.
- **Current three-cohort figures:** development evidence only until Feng and Zeller strict-production outputs replace historical inputs.

The final manuscript must report the exact completed sample and profile counts from the seals, not from plans, posters, Slurm job totals, or development manifests.

### 3.3 Critical work still required before submission

#### Priority 1 — required for valid final results

1. Complete and seal strict-production Feng and Zeller.
2. Audit exact profiler versions, databases, reference files, parameters, native abundance fields, and sample/profile counts.
3. Freeze the taxon identity and alias policy before rebuilding the final response atlas.
4. Build a single definitive three-cohort canonical input from the three sealed cohort products.
5. Run the frozen analysis plan once in `DEFINITIVE` mode.
6. Produce a participant/sample/profile flow table with explicit reasons for every exclusion or missing profile.
7. Regenerate every manuscript table and figure from the same definitive run.
8. Create a final source-data package connecting every plotted element to a machine-readable table.

#### Priority 2 — required for high scientific rigor

1. Report cohort-specific effects and uncertainty before pooled summaries.
2. Use biological-sample clustered or hierarchical uncertainty throughout repeated-dose analyses.
3. Report effect sizes and 95% confidence intervals, not only adjusted significance calls.
4. Make the multiple-testing family explicit for every endpoint.
5. Add sample-level bootstrap intervals for principal descriptive summaries where model-based intervals are not available.
6. Complete sensitivity analyses for zeros/pseudocounts, detection definition, native versus harmonized denominator, taxon aliases, sequencing depth, and reference/database representation.
7. Confirm that condition balance and sample completeness are preserved at every dose and target.
8. State whether community fractions are nominal or achieved and apply the frozen achieved-dose tolerance consistently.
9. Use a single frozen feature universe within each inferential family and document target exceptions.
10. Include a negative-control or falsification analysis. Strong candidates are phenotype-label permutation for biomarker fragility, shuffled target identities for response specificity, or spike-free null contrasts where technically available.

#### Priority 3 — analyses that would materially strengthen a Nature Communications submission

1. **Leave-one-cohort-out response transfer.** Learn target-to-feature response operators in two cohorts and predict responses in the unseen third cohort. This is a stronger generalization test than same-sample community reconstruction alone.
2. **Cross-cohort mixture reconstruction.** If component matching permits it, test whether independently learned response components improve community prediction in an unseen cohort relative to composition-only prediction.
3. **Truth-anchored biomarker fragility test.** Compare perturbation stability of known null discoveries generated by phenotype permutation with externally replicated and cohort-specific biomarkers.
4. **Stability-versus-transfer analysis.** Test prospectively whether perturbation-stable biomarkers replicate better across cohorts. If they do not, report that negative result: analytical stability and biological transportability may be distinct properties.
5. **Database/reference sensitivity.** Quantify how target genome representation and database version explain profiler-specific recovery, particularly for taxa with unexpected differences.

These additions should be prioritized by inferential value, not figure count. A coherent six-figure paper is stronger than a crowded manuscript containing every available analysis.

## 4. Recommended revised paper structure

### Proposed title

Preferred:

> **Controlled in silico perturbations reveal taxonomic profiling errors and biomarker fragility across colorectal cancer cohorts**

More general alternative:

> **Controlled metagenomic perturbations map profiler-specific measurement responses and downstream biomarker instability**

The title should avoid implying cellular ground truth. “Ground-truth” is defensible only when immediately qualified as read-level perturbation truth.

### One-sentence paper pitch

By adding known reads to 511 real stool metagenomes from three CRC cohorts, we show that taxonomic profilers have structured, taxon- and background-dependent response functions whose cross-talk and non-additivity propagate into biomarker stability and limit cross-cohort interpretation.

### Abstract structure

Use five compact moves:

1. **Problem:** profiler measurements influence microbiome biomarker discovery, but ordinary cohort comparisons lack feature-level analytical truth.
2. **Approach:** introduce matched read-level perturbations into real metagenomes across three cohorts, ten taxa, two perturbation designs, and two profiling workflows.
3. **Primary result:** detection and continuous recovery vary by profiler, taxon, dose, and clinical background.
4. **Mechanistic/downstream result:** perturbations generate structured non-implanted responses; single-target response information explains mixtures to different degrees across profilers; CRC biomarker calls show heterogeneous stability and limited cohort transfer.
5. **Meaning:** controlled perturbations provide a general measurement stress test and reliability annotation, but do not convert read fractions into biological abundance or directly label native biomarkers as true or false.

Do not put exact final numbers in the abstract until the definitive run is sealed.

### Introduction

Recommended four-paragraph logic:

1. CRC microbiome biomarkers are promising but replicate poorly across cohorts and pipelines.
2. Existing benchmarking commonly compares profilers without knowing the correct taxon-level response in complex real samples.
3. Read-level perturbations provide controlled analytical truth while retaining authentic community and technical backgrounds; clearly state what this truth does and does not represent.
4. State the revised questions:
   - How do profilers detect and quantify controlled taxon additions?
   - Do additions cause structured changes in non-implanted reported taxa?
   - Do independently measured responses predict community perturbations?
   - How do these measurement responses affect CRC biomarker stability and transfer?

The old filtering/calibration rationale should be shortened. It is now one possible application of the response framework, not the organizing thesis.

### Results

#### Result 1 — A three-cohort controlled perturbation framework preserves real sample backgrounds

Include:

- three cohorts and clinical composition;
- individual versus community designs;
- nominal and achieved fractions;
- two profilers and their native-output semantics;
- paired baseline design;
- exact final sample/profile flow;
- schematic plus compact quality-control panel.

Primary message: the design establishes read-level analytical truth in heterogeneous real metagenomes.

#### Result 2 — Detection and quantitative response depend on profiler, taxon, dose, and background

Primary endpoints:

- native detection probability;
- baseline-adjusted recovered signal;
- slope and deviation from the read-proportional expectation;
- continuous `log2(observed / expected)` as an interpretable secondary display.

Use recovery classes such as ≤10%, 10–50%, and >50% only as secondary descriptive summaries. The manuscript must not make a 10% cutoff appear biologically privileged.

The *F. nucleatum* case study can illustrate baseline detection and recovery, but it should not carry the general conclusion alone. Confirm its exact canonical identity and show that the result is not caused by suffix collapsing.

#### Result 3 — Controlled additions produce structured non-implanted taxonomic responses

Introduce the response signal:

`observed abundance - retained baseline abundance`.

For non-implanted taxa, estimate target-to-feature response slopes across dose. Demonstrate recurrence, direction, magnitude, and cohort dependence. Use “non-implanted taxa,” “collateral responses,” or “taxonomic cross-talk,” not “bystanders” if accessibility is a concern.

Interpret these as profiler-output behavior under a compositional sequencing intervention. Do not describe them as ecological interactions.

#### Result 4 — Independent response components explain community perturbations unevenly across profilers

Compare:

- a composition-only prediction, which retains the diluted baseline;
- an additive prediction formed from matched independent-spike response components.

The present poster analysis is a same-sample mechanistic reconstruction, not trained machine learning and not external validation. Say this explicitly. Report error ratios with uncertainty and stratify by cohort, dose, and whether the feature is an implanted target or non-implanted taxon.

If completed, add leave-one-cohort-out response-operator prediction as the stronger generalization test.

#### Result 5 — Controlled perturbations reveal heterogeneous stability of CRC biomarker inference

Fit the frozen CRC-versus-Control model to baseline and perturbed profiles using the same feature universe and sample panel. Separate direct implanted targets from non-implanted features. Report:

- retained baseline calls;
- lost calls;
- gained calls;
- direction reversals;
- continuous effect changes;
- ideal-counterfactual disease-dependent distortion.

Call these outcomes stability, fragility, or perturbation sensitivity. Do not call lost native biomarkers false positives.

#### Result 6 — Analytical stability and cross-cohort replication provide complementary evidence

Show the directional discovery-to-validation matrix across all three cohorts and both profilers. Rows are discovery cohorts and columns are independent validation cohorts. The diagonal reports discovery counts; off-diagonal cells report same-direction significant replication among evaluable biomarkers.

Then relate, without circularity, perturbation reliability to external replication. If stable biomarkers do not transfer better, make that a result rather than hiding it. It would support the conclusion that response stability and cohort transportability capture different failure modes.

#### Optional Result 7 — Response-informed calibration has limited, context-dependent transfer

Include this as a main result only if a frozen, training-only rule or calibration model improves held-out outcomes without compromising target recovery. The present development screen does not support a broad static blacklist. If the definitive analysis remains negative or asymmetric, place it in the supplement and use it to establish the limits of correction.

### Discussion

Recommended sequence:

1. State the conceptual advance: perturbation response is measurable and structured, and it affects downstream inference.
2. Explain why profiler differences arise: reference/database representation, classification strategy, marker dependence, compositional normalization, detection thresholds, and abundance estimation.
3. Discuss non-implanted responses and mixture reconstruction as evidence that analytical behavior is not captured by target detection alone.
4. Explain the biomarker implication carefully: robustness under controlled perturbation is a reliability dimension, not biological truth.
5. Discuss cross-cohort transfer and the distinction between analytical stability and population transportability.
6. State limitations prominently:
   - read-level rather than cellular or wet-lab ground truth;
   - only two profiling workflows and frozen database versions;
   - ten selected CRC-associated targets;
   - simulated reads and specific sequencing-error model;
   - no extraction, library preparation, or wet-lab effects;
   - dependence on available reference genomes and taxon aliases;
   - selected independent-spike subset, although frozen outcome-independently;
   - disease labels and covariates inherit limitations of the source cohorts.
7. End with the forward implication: controlled perturbation assays could accompany profiler/database releases and high-stakes biomarker studies as measurement validation, not as automatic biological correction.

### Methods

Suggested ordering:

1. Cohorts, metadata harmonization, eligibility, and ethics/source studies.
2. Frozen sample selection and analysis populations.
3. Read preprocessing and host removal.
4. Target selection, reference assemblies, and alias policy.
5. Read simulation and spike construction.
6. Achieved fractions and matched baseline definition.
7. Taxonomic profiling workflows, databases, and native abundance semantics.
8. Expected read-proportional reference and deterministic paired endpoints.
9. Detection model.
10. Continuous dose-response model.
11. Non-implanted response atlas.
12. Community superposition and prediction-error comparison.
13. Artificial differential-abundance model.
14. Native CRC biomarker model and frozen feature universe.
15. Ideal-counterfactual disease-distortion model.
16. Cross-cohort directional replication.
17. Assembly, denominator, zero, taxon-identity, and reference-representation sensitivities.
18. Multiplicity, uncertainty, missing-data policy, and software environment.
19. Reproducibility, code, data, manifests, and source-data availability.

Each method must state its experimental unit, estimand, sample panel, model formula, covariates, uncertainty estimator, multiplicity family, fallback policy, and whether it was prespecified or exploratory.

## 5. Proposed main-figure architecture

### Figure 1 — Study design and analytical estimands

Schematic of paired read perturbations in real metagenomes, two profiler workflows, individual/community designs, and the four output layers: detection, quantitative response, non-implanted response, and biomarker propagation. Add final cohort/profile counts.

### Figure 2 — Profiler- and taxon-dependent detection and quantitative recovery

Combine model-based detection curves or contrasts with continuous quantitative response. Keep the full species names and canonical order. The existing continuous community recovery plot is a useful panel, but should be accompanied by effect estimates and intervals.

### Figure 3 — Taxonomic cross-talk atlas

Show non-implanted response slopes by implanted target and profiler, with magnitude, direction, uncertainty, and recurrence across cohorts. A compact cohort-consistency panel can distinguish reproducible response structure from cohort-specific effects.

### Figure 4 — From individual responses to community perturbations

Show representative observed-versus-predicted examples and the error ratio comparing additive response-informed with composition-only prediction. Make the “better/equal/worse” direction immediately legible. Explicitly label same-sample reconstruction versus any held-out transfer panel.

### Figure 5 — CRC biomarker stability under controlled perturbation

Show target-excluded retention, loss, gain, sign reversal, and continuous effect change across dose. Add the ideal-counterfactual distortion result so that threshold crossings are not overinterpreted.

### Figure 6 — Cross-cohort biomarker transfer and its relationship to perturbation stability

Use the three-by-three directional replication matrix for both profilers. Add a panel testing whether perturbation reliability predicts external replication. If no association is observed, title the figure around complementarity rather than improvement.

### Supplementary figures

Place the following primarily in the supplement:

- *F. nucleatum* baseline prevalence and recovery case study;
- threshold recovery-class plots;
- minimum significant spike-fraction heatmap;
- full cohort/condition stratification;
- target-specific off-target patterns;
- profiler semantic and denominator audits;
- database/reference representation;
- assembly sensitivity;
- zero/pseudocount and prevalence-filter sensitivity;
- residual, calibration, convergence, and influence diagnostics;
- development/history comparison with the old unpaired analysis;
- calibration/filter Pareto curves and negative transfer findings;
- complete biomarker transition and response ledgers.

## 6. Statistical and interpretive standards

### Core estimands

For original native abundance `o`, total implanted fraction `F`, target fraction `f`, and observed post-implant abundance `a`:

- retained baseline: `(1 - F)o`;
- expected target abundance: `(1 - F)o + f`;
- expected non-target abundance: `(1 - F)o`;
- recovered target signal: `a - (1 - F)o`;
- response residual: `a - expected`.

These quantities should be introduced once and reused consistently.

### Required reporting behavior

- Present absolute denominators and sample counts with percentages.
- Report uncertainty around central estimates.
- Distinguish nominal from achieved dose.
- Distinguish absence from non-detection.
- Distinguish statistical non-significance from equivalence.
- Distinguish same-sample reconstruction from trained prediction.
- Distinguish development, exploratory, prespecified, and definitive analyses.
- Report negative and asymmetric transfer results.
- Preserve raw profiler identifiers in audit tables even when canonical names are shown in figures.
- Keep Control, Adenoma, CRC ordering everywhere.

### Language to use and avoid

| Prefer | Avoid unless directly justified |
|---|---|
| non-implanted taxon | bystander |
| collateral profiler response / taxonomic cross-talk | ecological interaction |
| perturbation-sensitive biomarker | false-positive biomarker |
| read-proportional expectation | true biological abundance |
| directional external replication | validation of biological truth |
| same-sample additive reconstruction | prediction model |
| native profiler abundance | directly comparable abundance between profilers |

## 7. Reproducibility and journal-readiness requirements

A high-tier submission should include:

1. immutable code release with DOI;
2. frozen container digest and complete software/database versions;
3. cohort manifests and sample inclusion/exclusion flow;
4. machine-readable analysis policy and model specifications;
5. checksummed source data for every main and supplementary figure;
6. statistical reporting checklist and computational reproducibility statement;
7. data-access instructions for controlled or externally hosted sequencing data;
8. a mapping from raw profiler identifiers to displayed canonical taxa;
9. seed, resampling, and permutation provenance;
10. an automated command or workflow that rebuilds the final manuscript outputs from sealed inputs.

Nature Communications evaluates whether conclusions are supported by data, appropriate controls, and rigorous analysis, as well as whether the work provides a substantial advance of broad relevance. Its current author resources also emphasize reporting summaries, data availability, code availability, and statistical/computational transparency. Before submission, check the current official [editorial process](https://www.nature.com/ncomms/submit/editorial-process), [submission guide](https://www.nature.com/ncomms/submit/guide-to-authors), [reporting resources](https://www.nature.com/ncomms/submit/resources), and [review criteria](https://www.nature.com/ncomms/for-reviewers/writing-your-report).

## 8. Claim-to-evidence gate

| Proposed claim | Minimum evidence required | Current status |
|---|---|---|
| quantitative recovery differs by profiler, taxon, dose, and background | definitive three-cohort paired models, CIs, sensitivity analyses | implemented; definitive rerun pending |
| controlled additions cause structured non-implanted responses | definitive response atlas, uncertainty, recurrence, identity sensitivity | implemented in development; definitive rebuild pending |
| independent responses improve community reconstruction | matched component audit, composition-only comparator, uncertainty | development result available |
| response behavior transfers across cohorts | leave-one-cohort-out prediction with locked model | not yet established as a central definitive result |
| native CRC biomarkers differ in perturbation stability | definitive frozen disease models and complete sample panels | implemented in development; definitive rerun pending |
| perturbation stability identifies more reproducible biomarkers | external replication or truth-anchored null test performed without leakage | not established; current result may be negative |
| spike-informed correction improves inference | held-out, training-only calibration with full model/BH refit and preserved target recovery | not established broadly |
| one profiler is universally superior | consistent evidence across all primary endpoints and sensitivity analyses | should not be claimed |

## 9. Practical revision sequence

### Phase A — finish evidence generation

1. Seal Feng and Zeller strict production.
2. Freeze aliases and profiler semantic audits.
3. Build the definitive three-cohort input.
4. Execute the frozen primary models and required sensitivities.
5. Run the leave-one-cohort-out analyses and truth-anchored fragility test if feasible.
6. Seal the final result bundle.

### Phase B — lock the scientific story

1. Write a one-page claims document containing only conclusions supported by definitive outputs.
2. Select six main figures according to the architecture above.
3. Move case studies, threshold summaries, and exhaustive diagnostics to the supplement.
4. Freeze terminology and visual encodings across all figures.

### Phase C — rewrite rather than patch

1. Rewrite the Results around response, cross-talk, mixtures, and biomarker propagation.
2. Rewrite the Introduction and Discussion to match that logic.
3. Rewrite the Methods from the frozen analysis documents.
4. Write the Abstract last, using final effect sizes and sample counts.
5. Conduct separate statistical, taxonomic-identity, reproducibility, and claims audits.

### Phase D — submission package

1. Deposit code and source data.
2. Complete journal reporting and data/code availability requirements.
3. Prepare a cover letter emphasizing the general measurement-science contribution, not only CRC-specific benchmarking.
4. Ask independent readers to challenge every causal, biological-truth, “false positive,” and generalization statement.

## 10. Bottom-line assessment

The revised project is potentially much stronger than the preprint. The key advance is no longer simply that Kraken2/Bracken and MetaPhlAn recover spike-ins differently. It is that controlled perturbations expose a structured measurement-response system and allow its consequences for mixtures and biomarker inference to be studied directly across real biological backgrounds.

The main obstacle is not a lack of analyses. It is the transition from extensive development work to one sealed, prespecified, internally consistent body of evidence. Completing that transition—and making appropriately limited claims—will matter more for a Nature Communications-level submission than adding further decorative figures or post hoc thresholds.
