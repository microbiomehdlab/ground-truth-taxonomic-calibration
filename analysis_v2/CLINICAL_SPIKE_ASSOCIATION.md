# Clinical spike association

Compare unchanged original controls with community-spiked adenomas, within each
cohort. This complements, rather than replaces, DA3's same-background experiments.
Use exactly DA1's eligible people, age/sex adjustment, native fraction abundance,
log2(1+a/1e-8) transformation, pinned MaAsLin2 1.18.0 LM and full clinical feature
family BH. Freeze the family and people across all doses; do not reselect them
using spike outcomes. Missing reported species are zero, not missing people.

There are 42 new fits: three cohorts × two profilers × seven total community
fractions (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1). The ten-member equal-weight
panel's nominal per-target fraction is total/10. Six existing unspiked DA1
comparisons are reused verbatim, not refitted. All doses are retained regardless
of significance. No new upstream profiling or simulated reads are generated.

Every fit retains full-family results and native MaAsLin2 models; an independent
OLS calculation audits each estimable group effect, standard error and p-value.
The report includes 420 spiked target results, 60 original target references,
nominal 95% effect intervals, raw p/full-family q, native prevalence and abundance
summaries. Intervals are not simultaneous or multiplicity-adjusted. Saved source
snapshots, container hashes and exact input profile hashes protect provenance.
Failed attempts remain for diagnosis; no failed fit enters a completed report.

This asks whether a known added signal becomes significant in a real clinical
comparison. It cannot establish that natural adenoma effects exist, prove that
profiling alone explains their absence, or define a universal abundance cutoff.
Community additions also change relative composition; non-target discoveries
are not automatically false positives. Retain the existing null/robustness
warnings when interpreting ordinary MaAsLin2 significance.

## Cluster submission

After these scripts are committed and pulled on the cluster:

```bash
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA1_REPORT="$PWD/work/da1_clinical_results_20261004T100340Z/REPORT"
export DA_INVENTORY_ROOT="$PWD/work/three_cohort_da_inventory_20261001T151231Z"
export CLINICAL_SPIKE_ROOT="$PWD/work/clinical_spike_association_$(date -u +%Y%m%dT%H%M%SZ)"
export CLINICAL_SPIKE_CONCURRENCY=30
bash analysis_v2/submit_clinical_spike.sh
```

Preparation, 42 parallel array tasks (up to 30 concurrently), and collection are
submitted together. No existing analysis outputs are modified. After completion,
inspect `REPORT/status.json` and verify `REPORT/SHA256SUMS`, then download REPORT.
Full-family results and native fits remain in `results/` on the cluster.
