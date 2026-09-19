# Ideal-counterfactual difference-in-differences

This prospective module tests whether a profiler departs from an ideal
read-proportional spike reference differently in CRC and Control samples.
It consumes the existing disease-model profile manifest and abundance-long
table plus the deterministic paired-endpoint table, which supplies the exact
total (`F`) and target-specific (`f`) implanted fractions.

For sample-level baseline abundance `o` and observed spiked abundance `a`, the
expected abundance is

- bystander taxon: `e = (1 - F)o`;
- implanted target: `e = (1 - F)o + f`.

This is a prespecified linear-response reference for profiler-native fractions,
not a claim that read fraction and reported taxonomic abundance are identical.
The fitted outcome is

`u = log2(a + 1e-8) - log2(e + 1e-8)`.

Within each cohort, study, analysis population, target, assembly arm, profiler,
and positive dose, the primary model is

`u ~ condition + scaled_age + sex`.

An invariant age or sex term is omitted and recorded rather than silently
estimated.

The reported coefficient is CRC versus Control. HC3 robust standard errors are
used, and BH correction is applied across the complete frozen species universe
within each exact model context. The universe follows the disease model:
species with at least 10% nonzero baseline prevalence plus every implanted
target. Baseline and perturbed profiles must contain the same biological
sample panel, and exact endpoint/profile identities and fractions are validated
before fitting.

The target and bystander roles are always reported separately. A bystander
coefficient whose corrected q value is at most 0.05 (and whose interval excludes
zero) is evidence of phenotype-dependent departure from the read-proportional
reference. It is a measurement-distortion
estimand, not proof that any baseline disease association is true or false.
Likewise, loss of significance in a separate disease model is
perturbation-sensitivity, not false-positive adjudication.

Run inside the frozen analysis environment:

```bash
Rscript analysis_v2/scripts/fit_ideal_counterfactual_did.R \
  --profile-manifest biomarker_profile_manifest.tsv \
  --abundance-long biomarker_abundance_long.tsv \
  --paired-endpoints paired_endpoints.tsv \
  --outdir ideal_counterfactual_did
```

The fail-closed project runner is `run_ideal_counterfactual_did.sh`. It takes
the sealed disease run so it reuses the exact abundance universe and sample
metadata, plus the paired endpoint table containing exact implanted fractions:

```bash
export DISEASE_RUN=/path/to/sealed/disease_run
export PAIRED_ENDPOINTS=/path/to/endpoints/paired_endpoints.tsv
export ANALYSIS_SIF=/path/to/frozen_analysis.sif
export ANALYSIS_STATUS=DEVELOPMENT_ONLY
export OUTDIR=/path/to/new/ideal_counterfactual_did
bash analysis_v2/run_ideal_counterfactual_did.sh
```

The output directory is required to be new or empty. Outputs include
feature-level HC3/BH results, sample- and model-panel audits, settings, session
information, checksums, and `SUCCESS`. Any missing pair, duplicate identity,
invalid fraction (`0 < f <= F < 1`), inconsistent profile path, invalid
covariate, or rank-deficient model causes a hard failure. Development runs may
record and use a complete common sample intersection. Definitive runs add
`--require-complete-panel` through the project runner and fail if any eligible
sample would be excluded. The runner also rejects development-marked source
inputs in definitive mode.

The feature-level model table is written by default. The much larger
sample-by-feature residual table is opt-in with `--write-sample-residuals` for
small audits and fixtures; omitting it avoids a potentially enormous redundant
intermediate in the full three-cohort run. Its row count remains recorded in
the summary even when the table is not materialized.

BH correction is within each target-by-dose model. Any claim about the presence
of distortion across the complete stress-test battery therefore requires a
separately frozen global or hierarchical decision rule; individual context q
values must not be mined to declare that “some distortion” exists.
