# Calibration-to-biomarker linkage

**Status:** implemented and fixture-tested before definitive execution.

This module links two previously sealed outputs without refitting either one:

1. paired quantitative endpoints describing response to implanted sequencing
   evidence; and
2. paired spiked-versus-matched-baseline artificial-biomarker calls.

The primary calibration quantity is the continuous read-perturbation response
ratio `R = recovered_spike_signal / implanted_target_fraction`. `R = 1` is
read-proportional response on a profiler's native abundance scale. It is not
cellular-abundance accuracy. Exact achieved fractions are retained, while joins
and summaries use the frozen nominal grid (`0.01%, 0.05%, 0.1%, 0.5%, 1%, 5%`).

For descriptive presentation, context-level median `R` is labelled
`under_response` below 0.8, `read_proportional_band` from 0.8 through 1.2, and
`over_response` above 1.2. The continuous result is primary; the 20% band is a
prespecified descriptive aid, not a biological equivalence margin.

One linkage row represents a cohort, study, population, phenotype background,
target, assembly arm, profiler, nominal dose, and q threshold. It contains
sample-level response summaries alongside target recall, target effect,
precision, and off-target burden. Descriptive association tables report
Spearman correlations; they are not causal estimates and do not replace the
paired dose-response or biomarker models.

Run with:

```bash
export ENDPOINTS_FILE=/path/to/sealed/endpoints/paired_endpoints.tsv
export ENDPOINTS_SUCCESS=/path/to/sealed/endpoints/SUCCESS
export PAIRED_RUN=/path/to/sealed/paired_biomarker_run
export ANALYSIS_SIF=/path/to/frozen_analysis.sif
export ANALYSIS_STATUS=DEVELOPMENT_ONLY  # or DEFINITIVE
export OUTDIR=/new/linkage/path
bash analysis_v2/run_calibration_biomarker_linkage.sh
```

Status inheritance is fail-closed. A development source cannot generate a
definitive linkage package. Outputs include joined and summarized tables,
figure-source TSVs, figures, diagnostics, captions, provenance, and checksums.
