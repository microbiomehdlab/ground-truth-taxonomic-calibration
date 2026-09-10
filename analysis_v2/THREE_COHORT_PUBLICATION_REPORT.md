# Three-cohort publication overview

This layer assembles presentation-ready overview figures from completed cohort
reports. It does not refit models. Cohort-specific estimates remain visible,
and independent and community perturbations remain separate facets.

Inputs are one or more complete cohort analysis packages with sealed
`reports/artificial`, `reports/disease`, and `reports/calibration_linkage`
subtrees. Development mode may combine temporary packages and always writes
`DEVELOPMENT_ONLY.txt`. Definitive mode requires exactly three non-development
packages, a sealed three-cohort synthesis, and the sealed Yachida assembly-
sensitivity report.

The primary outputs are complete combined TSVs, q <= 0.05 figure-source TSVs,
four PDF/PNG overview candidates, draft captions, diagnostics, a run manifest,
and `provenance/publication_report.sha256`. The overview covers artificial-
target recall, native disease-biomarker retention, calibration response ratio,
and off-target discovery burden. The detailed cohort reports and assembly-
sensitivity supplement remain authoritative for lower-level results.

Example development execution:

```bash
export COHORT_RUNS="/path/to/yachida_development:/path/to/feng_zeller_development"
export EXPECTED_COHORTS="yachida,feng,zeller"
export ANALYSIS_SIF=/path/to/ground_truth_analysis_v1.sif
export REPORT_STATUS=DEVELOPMENT_ONLY
export OUTDIR="$PWD/work/analysis_v2_three_cohort_report_dev_$(date +%Y%m%d_%H%M%S)"
bash analysis_v2/run_three_cohort_publication_report.sh
```

For definitive execution, set `COHORT_RUNS` to the three definitive package
roots, `REPORT_STATUS=DEFINITIVE`, `SYNTHESIS_RUN` to the sealed cross-cohort
synthesis, and `ASSEMBLY_REPORT` to the sealed Yachida assembly-sensitivity
report. The runner rejects development evidence and incomplete packages.
