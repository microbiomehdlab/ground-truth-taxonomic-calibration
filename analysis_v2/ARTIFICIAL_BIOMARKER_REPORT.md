# Artificial-biomarker recovery report

**Status:** implemented and fixture-tested before definitive execution.

This report addresses the direct controlled question: can differential-abundance
analysis identify the species whose reads were implanted? Each spiked library is
compared with the unmodified library from the same biological sample. Analyses
remain separate by cohort, phenotype background, target, assembly arm, profiler,
and dose. The target is therefore a prespecified artificial biomarker, while any
other significantly enriched species is counted as an off-target call.

The report consumes a sealed output from `run_paired_biomarker_propagation.sh`.
At the primary BH threshold of `q <= 0.05` it produces target recall, precision,
off-target burden, and paired target-effect figures. A labelled `q <= 0.10`
sensitivity analysis is retained in the source tables. It also reports the
minimum tested fraction at which each target is recovered. This is an empirical
minimum among the tested doses, not a continuous limit of detection.

Exact achieved fractions are preserved in the figure-source table. Summary
tables and graphics use the frozen nominal grid (`0.01%, 0.05%, 0.1%, 0.5%,
1%, 5%`) so integer read-allocation differences cannot split one experimental
dose into multiple groups. The minimum-dose table records both the nominal dose
and its exact achieved fraction.

Run with:

```bash
export PAIRED_RUN=/path/to/sealed/paired_run
export ANALYSIS_SIF=/path/to/frozen_analysis.sif
export REPORT_STATUS=DEVELOPMENT_ONLY  # or DEFINITIVE
export OUTDIR=/new/report/path
bash analysis_v2/run_artificial_biomarker_report.sh
```

A definitive report cannot inherit a development analysis. Every figure has a
TSV source; tables, diagnostics, captions, provenance, and checksums are sealed.
Interpret the target effect as recovery of controlled sequencing evidence, not
as cellular abundance, biomass, or extraction efficiency.
# Integrated primary and secondary reporting

`run_artificial_biomarker_report.sh` accepts the sealed pooled analysis in
`PAIRED_RUN` and, optionally, the phenotype-stratified analysis in
`SECONDARY_PAIRED_RUN`. Rows are explicitly labelled `pooled_primary` or
`phenotype_stratified_secondary`. If pooled rows are present, the main figures
show the pooled primary estimand; both scopes remain in the source tables.

The minimum-dose table reports both the first observed significant detection
and the first **sustained** detection (the lowest tested dose for which the
target is detected at that dose and every higher dose). The former describes
sensitivity; the latter is less vulnerable to an isolated non-monotone call.

Where the model call table contains the full context columns, the report also
writes `off_target_call_ledger.tsv` and `recurrent_off_target_taxa.tsv`. These
preserve the identities of enriched non-target taxa, rather than reporting
only their count.
