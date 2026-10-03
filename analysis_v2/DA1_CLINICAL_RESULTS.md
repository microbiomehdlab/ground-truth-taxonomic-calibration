# DA1: actual clinical associations

This export completes the computational clinical-results table, not a claim
that any target must be associated with disease. It reuses the 12 existing
MaAsLin2 1.18.0 fits from the clinical pilot: three cohorts × two profilers ×
two contrasts (Adenoma versus Control and CRC versus Control). Those contexts
already used all eligible clinical samples and the complete frozen families.
The original files and their engineering status are preserved unchanged.

## Fixed method and interpretation

- Native profiler abundance fractions; `log2(1 + abundance / 1e-8)`.
- Unpaired clinical LM: disease group + age in years + sex; Control=0,
  disease=1, sex reference Female. Complete-case age/sex eligibility.
- MaAsLin2 normalization NONE, transformation NONE (input already transformed),
  standardize FALSE; no additional backend prevalence/variance filter.
- Baseline prevalence family frozen in the checkpoint, union ten targets.
  Two-sided group p-values; BH separately within each full context family,
  including non-estimable bookkeeping p=1. Raw p remains NA for non-estimable
  features. Native package q-values are retained separately, not substituted
  for full-family group-only BH.
- Positive and negative associations both retained. The existing
  `positive_discovery` field is directional, not the complete association set.
- `beta` is an age/sex-adjusted difference on this transformed scale, **not a
  raw-abundance fold change**. Confidence intervals are nominal 95% ordinary-LM
  t intervals (n−4 residual degrees of freedom), not multiplicity-adjusted or
  heteroscedasticity-robust intervals. Previously identified model limitations
  and sensitivity analyses remain relevant; no new method is selected here.
- Prevalence is the fraction with positive reported native abundance among
  eligible people, not proof of biological presence/absence. Distributions
  are descriptive and unadjusted, whereas beta is age/sex-adjusted.
- A nonsignificant association does not demonstrate biological absence or
  prove a profiling cause. Interpret alongside spike recovery and DA3.

## Verification and outputs

The exporter verifies the sealed pilot plan and every original fit-bundle
checksum, exact context/sample/family identities, age/sex eligibility, the
pinned package session, native group statistics, estimability and full-family
BH. Independent OLS arithmetic checks the saved fit against its specified
design; the published beta/SE/p/q remain the **saved MaAsLin2 values**, not
replacement inference. Intervals and eligible-group distributions are added.
There is no automatic fallback or refitting. An inconsistency retains an
attempt directory with FAILED.json and prevents final report publication.

`REPORT/` contains:

- `target_results.tsv`: 120 target–clinical-context rows, including failures
  of estimability and nonsignificant results.
- `full_family_results.tsv`: every frozen feature, not only significant taxa.
- `context_summary.tsv`: 12 contexts, group sizes, exclusions, family sizes,
  estimability, significant counts and residual degrees of freedom.
- `eligibility.tsv`: all checkpoint candidate people and exclusion reasons;
  controls appear in both contrasts, so rows are not unique-person counts.
- `contexts/`: reproducible saved abundance/metadata/context/results and source
  provenance. `original_bundle_checksums.txt` describes the full original
  bundle; it is NOT a checksum manifest for this reduced copied subset.
  Large native models remain in the original verified pilot directories.
- `source_evidence.tsv`, `status.json`, `SUCCESS`, `SHA256SUMS`.

Source code and container identity are retained in the enclosing execution
directory. Report PASS means successful audited export. Scientific review and
final figure/manuscript interpretation follow; DA2/DA3 are separate workflows.

## Cluster execution

After the implementation is committed and pushed, pull the branch and run:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA_PILOT_PLAN="$PWD/work/da_pilot_plan_20261001T161335Z"
export DA_PILOT_RESULTS="$PWD/work/da_pilot_results_20261001T161641Z"
export DA1_EXPORT_ROOT="$PWD/work/da1_clinical_results_$(date -u +%Y%m%dT%H%M%SZ)"
bash analysis_v2/submit_da1_clinical_export.sh
```

One 1-CPU/16-GB job audits all contexts. It snapshots committed code and only
reads existing data; it does not affect running DA3 jobs. On completion read
`$DA1_EXPORT_ROOT/REPORT/status.json` and verify its `SHA256SUMS`. Keep original
fit bundles until the public release retention review; do not delete models.
