# Reporting layer

The reporting layer converts sealed analysis tables into manuscript-source
tables, figure-source data, figures, diagnostics, draft captions, provenance,
and checksums. Figures are views of preserved TSV sources; they are not the sole
record of any result.

`run_disease_biomarker_report.sh` currently implements the disease-biomarker
module. It produces target recall, biomarker-set stability, and target disease-
effect-change displays at the primary q <= 0.05 threshold, while preserving
both prespecified thresholds in source tables. It also exports a model-covariate
audit, which makes invariant-covariate omission visible.

Status inheritance is fail-closed: a source containing `DEVELOPMENT_ONLY.txt`
cannot produce a `DEFINITIVE` report. Development figures are for integration
and presentation testing only. The definitive package will be regenerated from
sealed complete cohort analyses without changing code or graphical definitions.
