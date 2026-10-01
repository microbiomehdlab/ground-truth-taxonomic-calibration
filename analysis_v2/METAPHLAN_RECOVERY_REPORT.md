# MetaPhlAn descriptive recovery reference report

Input: verified three_cohort_metaphlan_recovery_20261001T140855Z,
cluster job 3097831, PASS in 36 seconds, empty stderr; root checksums verified
by André. Select only profiler=metaphlan4: expected 41,170 positive-dose rows.
Bracken rows in that bundle are not used by this report.

summarize_metaphlan_recovery.py compares response_ratio_profiler_scale
(genome-equivalent primary) against response_ratio (read-fraction sensitivity)
within the same physical target endpoint. Raw abundances are unchanged.
Groups: cohort, population, clinical condition, target, nominal total dose.
Nominal dose is decoded from the canonical builder's profile ID _f tag;
achieved total/target fractions remain separate range columns.

Summaries include median, Q1, Q3, minimum, maximum, negative/zero recovery
counts and native non-detection counts. Absolute relative error is |ratio-1|.
Paired deltas are calculated per sample before summarization, not differences
of independently computed medians. Preserve all zero and negative responses.
IQRs are descriptive sample variation, not confidence intervals. A smaller
error under the revised reference is not evidence of improved profiler output.

run_metaphlan_recovery_report.sbatch uses the shared Bracken report plotting
script with explicit MetaPhlAn labels selected from summary.json. Existing
Bracken behavior remains its default. Outputs: two summary TSVs, source hashes,
summary.json, 18 PDF and 18 PNG figures, R session info, root code/image/commit
provenance, SUCCESS and SHA256SUMS. Tables and figures use distinct directories.

Required exports: PROJECT, ANALYSIS_SIF, METAPHLAN_RECOVERY_ROOT and fresh
MPA_REPORT_ROOT. Tests: test_metaphlan_recovery_report.py (3) and
test_bracken_recovery_report.py (6); local synthetic rendering passed.
The synthetic minimal plot fixture has one dose and one target, so line-group
warnings and empty facets are expected there, not evidence about real data.
Run tests inside the pinned image on the cluster. Validate completed report
checksums before copying it locally for scientific interpretation.
