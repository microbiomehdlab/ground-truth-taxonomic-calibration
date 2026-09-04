# Profiler units and denominator decisions

**Status:** command audit complete; representative final-output audit pending.

## Frozen native outputs

Kraken2 classifies paired FASTQs and Bracken estimates species counts with
`-l S`, read length 100 in production, and abundance threshold 10. Its native
species file contains `new_est_reads` and `fraction_total_reads`. The v2
analysis retains both fields and does not silently close the species table to
100%. A remainder in the species fractions is not automatically labelled
“unclassified”: it can also reflect ranks or taxa absent from the retained
Bracken species table and the configured threshold.

MetaPhlAn 4 profiles paired FASTQs with the frozen database and
`--ignore_eukaryotes --ignore_archaea`. The command does not request
`--unclassified_estimation`. Its native `relative_abundance` percentage is a
marker-based compositional estimate; it is not an input-read assignment
fraction. Species rows, any explicit `UNCLASSIFIED` row, and non-species rows
are audited separately.

## Decisions for v2

- Native profiler abundance is the primary representation.
- Baseline-adjusted change and dose-response are computed within profiler.
- No primary calculation converts MetaPhlAn percentage into read counts.
- No primary calculation assumes the complement of species rows is exactly
  unclassified for either profiler.
- Non-detection remains zero in unconditional summaries. A detected-only
  analysis, if used, is secondary and clearly conditional.
- Species-closed renormalization may be reported only as a labelled sensitivity
  analysis; it is not “denominator harmonization.”
- There is currently no single abundance denominator that makes the two native
  outputs biologically identical. Cross-profiler conclusions concern response
  to the same implanted read perturbation, not equality of measurement units.

## Legacy caveats (do not rewrite historical results silently)

`workflows/kraken2_bracken/postprocess_local.sh` copies its closed-world CSV to
the historical open-world filename. It therefore does not create an independent
open-world representation. `workflows/metaphlan4/postprocess_local.sh` creates
pseudo-read counts by multiplying relative-abundance percentages by FASTQ read
counts and labels the complement of species percentages as unclassified. Those
derived values must not be treated as native assigned reads or as a measured
unclassified fraction in v2.

The paired QC helpers also count R1 plus R2 FASTQ records. Before any comparison
to Kraken2/Bracken fragment counts, read-versus-pair units must be verified; a
factor-of-two unit mismatch must not be allowed into a harmonized endpoint.

## Required empirical gate

Run `scripts/audit_profiler_semantics.py` on representative baseline and spiked
native files from every cohort, including a detection and non-detection case.
Archive its TSV, input checksums, and `SUCCESS` marker with the downstream run.
Freeze detection and transformation rules only after this audit, without using
the direction of final profiler differences to choose them.
