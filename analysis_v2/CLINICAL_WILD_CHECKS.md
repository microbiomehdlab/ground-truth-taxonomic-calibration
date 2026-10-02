# Clinical wild-bootstrap engineering candidate

Separate from DA3 and all primary results. No automatic production authorization.
The clinical group coefficient remains OLS, adjusted for age and sex. The candidate
uses a group-null restricted fit, HC2-adjusted restricted residuals, independent
person-level Rademacher multipliers, and HC3-studentized bootstrap statistics.
Common person weights across features retain residual dependence. This is an
approximate conditional bootstrap, not an exact randomization test or proof of
clinical FDR control. The plus-one Monte Carlo p-value uses 9,999 draws.

All twelve cohort/profiler/clinical-contrast contexts are covered. The actual-data
comparison uses the first 20 lexically sorted frozen features, selected without
looking at effect sizes or significance. It is an engineering subset, not a
representative target panel and not the manuscript's discovery family. Gaussian,
skewed and group-heteroskedastic global-null checks use the same actual age/sex
design and 20 synthetic features, 200 independent simulations per context and
scenario. Each method receives identical simulated outcomes. Each scenario's 200
simulations are split into four independent jobs. There are 156 tasks total.

The minimum Monte Carlo p is 0.0001, below 0.05/20. It is NOT below 0.05/3471:
full-family production needs a separately resolved precision and cost policy.
Null rejection rates are not biological results. Fewer discoveries do not itself
establish better validity; inspect pointwise rejection, family rejection, uncertainty,
positive controls and full-family feasibility before choosing a primary method.

The launcher snapshots committed code, records image checksum and job IDs, and
uses at most 30 concurrent one-CPU tasks by default. Workers verify image and
input checksums. Failed tasks remain failures, not zero-discovery results. The
after-any collector reports missing/failed tasks and exits nonzero if incomplete.
It never changes source pilot results or the running DA3 snapshot.

Local checks:

```bash
Rscript analysis_v2/tests/test_clinical_wild_bootstrap.R
Rscript analysis_v2/tests/test_clinical_hc3.R
python3 -B analysis_v2/tests/test_clinical_wild_checks.py
bash -n analysis_v2/submit_clinical_wild_checks.sh
bash -n analysis_v2/run_clinical_wild_checks.sbatch
```

Cluster submission (after pulling the commit containing these files):

```bash
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA_PILOT_RESULTS="$PWD/work/da_pilot_results_20261001T161641Z"
export WILD_ROOT="$PWD/work/clinical_wild_checks_$(date -u +%Y%m%dT%H%M%SZ)"
export WILD_CONCURRENCY=30
bash analysis_v2/submit_clinical_wild_checks.sh
```

Read `jobs.tsv` to recover IDs after reconnecting. When the report job completes,
inspect `REPORT/status.json`, `REPORT/summary.tsv` and task `model.err` files.
Checksummed reports and per-task outputs are retained for review.
