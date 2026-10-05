# Prespecified conversion-factor stress test

Decision recorded before examining new model outcomes: use multipliers 0.5 and
2.0 on the baseline G_eff entering MetaPhlAn's reference conversion. Retain the
completed primary multiplier 1.0 results without refitting. These deliberately
broad symmetric-on-log-scale scenarios are stress tests, NOT plausible ranges
estimated from genome-size uncertainty, confidence intervals, or alternative
profiler measurements. They uniformly scale G_eff/G_target; they do not test
taxon-specific genome-size errors, marker breadth, or baseline measurement bias.

Scope: MetaPhlAn only, the same 360 original pilot contexts, three cohorts, n=10/20,
three tested doses, uniform exposure, first 20 saved allocations. Two scenarios
give 720 new reference contexts, 7200 new target/scenario rows, 72 ten-context
Slurm tasks. Combined report includes primary + two stress scenarios: 10800 rows.
No new profiling, new biological samples, new observed MaAsLin2 fits or change to
the declared transform, filtering, native scale, families or full-family BH.

For each added library: q_t=F*(N_t/N_total)*(multiplier*G_eff/G_t),
D=(1-F)+sum(q_t). Whole-profile background dilution and target addition are
recomputed with that D. Zero-dose reference profiles remain the original baseline.
Use verified integer counts and exact achieved fractions from reference_audit.tsv.
Native inputs and observed primary results are unchanged. The original snapshotted
worker, R fitting code and analysis image are reused, with checked image identity.
New orchestration/collection is snapshotted independently. Existing outputs are
never edited; preparation and completed-task collection enforce checksum identity.

## Interpret all outcomes

Compare coefficients, SE, raw p, full-family q, primary/HC3 categories and NE, not
only significance. Original observed all-zero/perfect-fit classifications cannot
change through reference scaling. What can change is how often the stipulated
reference is statistically significant. A persistent gap strengthens robustness
to this one conversion assumption, not causal identification of natural disease
failure. A disappearing gap limits the expected-detectability claim; retain it.
Do not pick a multiplier to maximize discoveries or change the reported primary.

## Cluster submission

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export REFERENCE_ROOT="$PWD/work/reference_response_pilot_20261005T202107Z"
export REFERENCE_SCALE_ROOT="$PWD/work/reference_scaling_$(date -u +%Y%m%dT%H%M%SZ)"
export REFERENCE_CONCURRENCY=72
bash analysis_v2/submit_reference_scaling.sh
```

Three submissions: preparation, array of 72 tasks, dependent collector (74 tasks
total). Concurrency controls tasks, not dedicated nodes; each worker requests
1 CPU / 12 GB. Choose an allowed cluster resource budget. Empty cluster is not
permission to exceed scheduler limits. Preparation failure blocks dependents;
inspect the preparation error log, do not remove dependencies.

```bash
cat "$REFERENCE_SCALE_ROOT.jobs.tsv"
cat "$REFERENCE_SCALE_ROOT/REPORT/status.json"
(cd "$REFERENCE_SCALE_ROOT/REPORT" && sha256sum -c --quiet SHA256SUMS)
```

Download only REPORT for scientific review. Full-family variant results remain
on cluster. No automatic production authorization or further scope expansion.

## Verification scope

Local sealed synthetic fixture covers preparation, all 360-context scope/counts,
half/double construction, unchanged zero-dose baseline and observed targets,
snapshotted worker routing, 10800-row combined collection, source-tamper refusal
and fresh-output guards. Existing reference formula and NE fixtures also pass.
Real MaAsLin2 backend was validated in the original cluster pilot, but the new
sensitivity has not been executed on real data locally. Original workers retain
full-family coverage checks, numerical equivalence checks and atomic resume.
