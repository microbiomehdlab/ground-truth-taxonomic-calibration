# Cross-cohort disease-biomarker synthesis

Yachida, Feng, and Zeller remain separate studies. Cohort-specific disease
estimates are the primary evidence; pooling is a summary and never substitutes
for displaying those estimates.

For each shared population, target, assembly arm, profiler, nominal dose,
disease contrast, feature, and model specification, the synthesis uses a REML
random-effects model. Uncertainty uses modified Hartung–Knapp inference with a
scale bounded below by one. The output includes tau-squared, Cochran's Q, I²,
a 95% prediction interval, direction consistency, and leave-one-cohort-out
estimates. With only three cohorts, heterogeneity estimates and prediction
intervals are inherently imprecise and must be interpreted cautiously.

Only a feature with a finite positive standard error in every prespecified
cohort is pooled. `synthesis_coverage.tsv` records incomplete or ineligible
cells. Dose alignment uses the prespecified dose label; the median and range of
realized target fractions are retained. BH correction is across features within
each population/target/arm/profiler/dose/contrast/model synthesis family.

The script is `scripts/meta_analyze_disease_biomarkers.R`. Definitive execution
requires all three sealed cohort-specific primary result tables. Development is
limited to synthetic fixtures until those tables exist.
