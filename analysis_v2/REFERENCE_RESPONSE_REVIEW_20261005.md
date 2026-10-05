# Reference-response pilot: paper handoff, 5 October 2026

Update: the targeted NE audit is now executed and locally verified. See
REFERENCE_NE_REVIEW_20261005.md: MetaPhlAn constants are all-zero inputs; Bracken
perfect fits are zero controls versus identical positive exposed values. The
first next-step input inspection below is completed, not still pending.

Reference assumptions and compact-figure scope are now reviewed in
REFERENCE_ASSUMPTIONS_REVIEW_20261005.md. Four-panel draft and source tables live
outside Git in ../reference_response_compact_20261005_v2/. No expanded reference
jobs were submitted; the next methodological step is construction sensitivity,
not another statistical-method search or a rerun of completed NE audits.

Status: executed and checksum-verified; initial descriptive scientific review
completed. Figure drafts exist. Not a production authorization or final paper
claim. This entry supersedes the earlier statement that this pilot still awaited
cluster execution; other pending analyses are not marked complete.

## Evidence and methods

Cluster run: `work/reference_response_pilot_20261005T202107Z`.
Preparation 3108874 and all 72 array tasks 3108875 completed. Collector 3108876
failed because the output overlap guard rejected its designated REPORT child.
Report-only fix f4b8edb preserved sealed calculations; retry 3108948 completed
in 24 seconds, exit 0:0. REPORT/status.json reports 720 contexts and 7200 target
rows, COMPLETE_PENDING_SCIENTIFIC_REVIEW, production_authorized=false. Package
checksums were verified on cluster and locally.

Three cohorts × two profilers × n=10/20 per artificial group × three nominal
total community doses (0.0001/0.001/0.01) × 20 saved allocations = 720 contexts;
ten targets/context. Adenoma backgrounds, uniform additions, distinct artificial
case/control people. Nominal per-target percentages are 0.001/0.01/0.1%.
These comparisons are not adenoma-versus-control clinical associations.

Pinned MaAsLin2 1.18.0 primary, positive effect and full-family BH q≤0.05,
log2(1+native abundance/1e-8), HC3 sensitivity. Baseline-anchored expected
response is not unbiased community truth. Construction and scale limitations:
REFERENCE_RESPONSE_EXTENSION.md. No new fits were required for the local review.

Local evidence is outside Git in sibling directories:

- `../reference_response_pilot_20261005T202107Z_REPORT/`: original report.
- `../reference_response_review_20261005/`: 30-page PDF, 30 PDF/SVG/450-dpi PNG
  panels, accessibility review PDFs, summaries, notes and checksums.

Original report SHA256SUMS digest:
`a35200c7dd888f39272aa12856ac062c2001c4b87c0c73484116de82cd7bbdd1`.
Scripts: summarize_reference_response.py, plot_reference_response.R;
definitions/regeneration: REFERENCE_RESPONSE_FIGURES.md. Private manuscript's
shared theme was used, without committing the LaTeX project. Representative
significance and effect pages visually inspected; exhaustive visual acceptance
remains pending. IQRs across allocations are not confidence intervals.

## Initial findings, including caveats

At nominal 0.001% per target and n=20/group, pooled 600 target-allocation rows/tool:
Bracken observed positive 360, reference positive 457; MetaPhlAn observed 1,
reference 495. Each arm's own flag is counted, including when the other is NE.
Do not infer these totals solely from BOTH_POSITIVE/REFERENCE_ONLY categories.
The low-dose observed/reference gap is particularly marked for MetaPhlAn;
higher doses show more agreement. Full target-level tables retain heterogeneity.

Bracken's Bfrag and Csym both have 0/60 reference-positive comparisons at that
setting. Some expected low-dose responses remain insufficient under this design;
profiling therefore cannot explain all absent significance.

All 148 observed Bracken NE rows are NON_ESTIMABLE_PERFECT_FIXED_FIT, not evidence
of missed spikes. All 291 observed MetaPhlAn NE rows are NON_ESTIMABLE_CONSTANT;
the report alone does not establish whether the constant is zero. No reference
NE rows occur. Effects/SE summaries use only jointly estimable allocations and
print their count. Never merge these reasons into one measurement-failure rate.
Primary versus HC3 category differs in 128/7200 comparisons (1.8%); this is not
a certificate of randomization validity, calibration or natural disease power.

## Next steps, in order

1. Inspect the saved inputs for constant MetaPhlAn cases and perfect-fit Bracken
   cases: zero/nonzero values, transformed vectors, and native/reference changes.
   Use existing inputs/results, not a new DA-method search or reprofile.
2. Review reference construction sensitivity, effect/SE results and each figure's
   denominators; decide whether the pilot supports expansion to all 100 saved
   allocations and U/P25/P50/P75. No scope expansion based on favorable results.
3. Choose a compact observed/reference + NE main figure, retaining all taxa/cohorts
   in supplements. Move crowded annotations if needed before figure acceptance.
4. Draft a Methods subsection and cautious Results subsection in the private
   manuscript, linking these outputs to recovery, baseline heterogeneity and
   clinical-reference evidence. Do not claim causal explanation of absent natural
   adenoma biomarkers or call artificial groups disease-association estimates.
5. Commit reviewed scripts/documentation separately from private manuscript/data;
   later freeze the final release and list the selected source tables for Zenodo.

No new cluster run is required simply to discuss these completed pilot plots.

### Targeted NE follow-up prepared

The user requested the next review step. A saved-input audit is implemented in
audit_reference_nonestimable.py, with sbatch entry point and four passing local
fixtures; see REFERENCE_NONESTIMABLE_AUDIT.md. Real input caches were not included
in the downloaded report, so actual zero/nonzero conclusions await this read-only
cluster audit. It performs no profiling/model refits and preserves original seals.
