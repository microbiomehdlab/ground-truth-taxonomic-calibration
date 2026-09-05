# Native disease-biomarker model

This model answers the manuscript-facing question that the paired perturbation
model cannot answer: how controlled read implantation changes actual
CRC-versus-Control and Adenoma-versus-Control biomarker conclusions.

Models are fitted separately by cohort, study, analysis population, target,
assembly arm, profiler, and dose. The biological sample is the unit of
inference. Each family must contain exactly the same samples at baseline and at
every positive dose, with Control, Adenoma, and CRC all represented.

For every native species abundance, the primary model is
`log2(abundance_fraction + 1e-8) ~ condition + scaled_age + sex`. The primary
contrast is CRC versus Control; Adenoma versus Control is secondary. HC3 robust
standard errors protect against heteroskedasticity. A complete-case model adding
scaled BMI is a prespecified sensitivity analysis because BMI missingness is not
uniform across cohorts. A covariate that is invariant in an analysis population
is non-estimable and cannot confound its within-population contrast; it is
omitted automatically and recorded in every result row. Other rank deficiency
remains a hard model failure rather than triggering silent term selection.

The species universe is fixed across doses within a family. Species require at
least 10% nonzero prevalence over all profiles in that family; the implanted
target is retained regardless of prevalence. BH correction is performed across
tested species within each cohort/study/population/target/arm/profiler/dose/
contrast/model family. The baseline call set is an observed disease contrast.
Disease propagation is therefore summarized by retained, lost, and gained
biomarkers, baseline retention, Jaccard stability, direction flips, and effect
changes. Because targets are implanted across phenotype groups, target
significance is a spurious-association diagnostic and is never called recall.
Jaccard stability and baseline retention are undefined when no baseline or
union call set exists; empty sets are not presented as perfect stability.

Profiler-native outputs remain distinct practical estimands. The model does not
claim that Bracken and MetaPhlAn measure identical cellular abundance. It tests
how each profiler's reported species abundance propagates into disease-marker
conclusions under the same controlled sequencing perturbation.

The fail-closed entry point is `run_disease_biomarker_propagation.sh`. It
requires a validated canonical table, a frozen sample metadata manifest, the
frozen analysis image, a new output directory, and an explicit
`DEVELOPMENT_ONLY` or `DEFINITIVE` status. Development results must not be used
as manuscript evidence.
