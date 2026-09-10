# Reporting layer

The reporting layer converts sealed analysis tables into manuscript-source
tables, figure-source data, figures, diagnostics, draft captions, provenance,
and checksums. Figures are views of preserved TSV sources; they are not the sole
record of any result.

`run_disease_biomarker_report.sh` currently implements the disease-biomarker
module. It produces baseline-marker retention, biomarker-set stability, and target disease-
effect-change displays at the primary q <= 0.05 threshold, while preserving
both prespecified thresholds in source tables. It also exports a model-covariate
audit, which makes invariant-covariate omission visible.

`run_artificial_biomarker_report.sh` reports direct recovery of the implanted
target. `run_calibration_biomarker_linkage.sh` links those calls to quantitative
read-perturbation response, preserving exact achieved fractions and producing
separate figure-source, diagnostic, caption, and provenance files.

`run_combined_calibration_report.sh` consolidates the cohort-specific detection
and continuous-recovery fits across Feng and Zeller, retaining independent and
community spike experiments as separate facets. It writes the complete model
contrasts, figure-source TSVs, combined diagnostics, draft captions, and sealed
PDF/PNG candidate figures. The same entry point is used for development and
definitive inputs; status inheritance prevents development results from being
presented as final evidence. Low doses are displayed on an ordered categorical
axis, detection uncertainty uses interval bars, and interpretation tables state
the proportional-recovery classification and cross-cohort directional
consistency without replacing the underlying estimates.

Status inheritance is fail-closed: a source containing `DEVELOPMENT_ONLY.txt`
cannot produce a `DEFINITIVE` report. Development figures are for integration
and presentation testing only. The definitive package will be regenerated from
sealed complete cohort analyses without changing code or graphical definitions.

`run_three_cohort_publication_report.sh` is the final overview layer. It
combines sealed cohort reports without refitting models, preserves cohort and
independent/community facets, and exports artificial-target recall, native
disease-marker retention, quantitative response, and off-target burden panels
with their complete source TSVs. Definitive mode requires exactly three
non-development cohort packages plus sealed cross-cohort synthesis and Yachida
assembly-sensitivity reports. `MANUSCRIPT_RESULTS_CHECKLIST.tsv` maps each
planned claim to its estimand, source evidence, candidate display, and gate.
The input contract and execution examples are in
`THREE_COHORT_PUBLICATION_REPORT.md`.
