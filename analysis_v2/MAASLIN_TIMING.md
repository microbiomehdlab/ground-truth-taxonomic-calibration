# Direct MaAsLin2 timing trial

Requested 3 October 2026. Purpose: measure whether calling the actual pinned
MaAsLin2 package for the existing DA3 design is affordable. This is not a new
DA method, a clinical association analysis, or automatic authorization of a
120,000-context rerun. Existing results and inputs are read-only.

## Fixed design and settings

48 contexts: three cohorts × Adenoma/CRC backgrounds × two profilers ×
n=5/20 per independent group × uniform community doses 0.0001/0.001.
These are total mixture read fractions, not per-taxon doses. Select the first
context ID in each design stratum from the existing sealed 120,000-context
plan, without reading its results. Use its exact participants, profiles and
entire frozen feature family. This deliberately measures small and larger
experiments and weak/moderate additions; it is not a representative random
sample of all scenarios or allocations.

Real backend: **MaAsLin2 1.18.0**, group-only unpaired LM, one core,
externally transformed `log2(1 + native_fraction/1e-8)`, normalization NONE,
transform NONE, standardize FALSE, no additional abundance/prevalence/variance
filter. Keep native tables and saved models; disable only plot generation.
These are our existing configured settings, **not MaAsLin2 defaults**. This
does not change the age/sex-adjusted clinical DA1 model or paired DA2 analysis.

The existing wrapper accounts for constants/perfect fixed fits as
non-estimable (raw p=NA, multiplicity bookkeeping p=1). Compare coefficients,
estimability, p-values, full-frozen-family BH q-values and positive discoveries
against the existing fast group-only LM, on exactly the same native matrix.
Numeric agreement uses absolute 1e-8 plus relative 1e-6 tolerance; significance
calls must match exactly. Native package q-values are retained separately:
they can differ from full-family q when features are non-estimable.
Warnings or unexpected missing fits fail a task and retain diagnostic files.
Agreement validates implementation, not type-I-error calibration; previous
null and heteroscedasticity limitations are not erased.

## Execution

The launcher snapshots committed code and records the container hash. It
submits preparation (including a mandatory real-backend fixture), a 48-task
array, and automatic collection. Dependencies use afterany so missing tasks
are reported, not silently excluded. No tasks submit a larger run.

For André to run on the cluster:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA3_BATCH_PLAN="$PWD/work/preproduction_checks_20261002T192903Z/draft_da3_plan"
export MAASLIN_TIMING_ROOT="$PWD/work/maaslin_timing_$(date -u +%Y%m%dT%H%M%SZ)"
export MAASLIN_TIMING_CONCURRENCY=12
bash analysis_v2/submit_maaslin_timing.sh
```

After collection, inspect `REPORT/status.json`, `cost_by_profiler_and_n.tsv`
and `rough_projection.tsv`; verify `REPORT/SHA256SUMS`. `jobs.tsv` preserves
job IDs across shell sessions. A report collection SUCCESS does not override
failures or mismatches in status.json. Failed attempts remain untouched; choose
a new root for a rerun. Download the small REPORT first, retaining per-task
native models on the cluster for any necessary mismatch inspection.

## Cost interpretation and validation status

`maaslin_seconds` includes package loading, wrapper checks, fitting and native
output I/O. `fast_seconds` averages 50 in-memory calls including transformation;
it excludes R startup and profile I/O. `task_seconds` additionally includes
profile loading/checking and R startup, but not Slurm wait, container startup
or final sealing. Threads are fixed to one. Runtime is measured at the chosen
array concurrency, so shared-storage contention may affect it.

The full-grid projection interpolates n=10/15 between n=5/20, averages the two
tested doses within each cohort/background/profiler, and weights by the actual
plan counts. It assumes untested arms/allocations have similar costs. Ideal
30-worker wall time is arithmetic, not a scheduler prediction or guaranteed
deadline. Batching may reduce startup cost; retaining every native model can
make disk usage substantial. Decide after reviewing time, disk and equivalence.

Local tests cover deterministic selection, immutable input checks, identity
and overlap rejection, failed-backend preservation, comparison and collection.
Local MaAsLin2 is unavailable: actual package execution is mandatory in the
cluster preparation job; real-input timing has not yet been measured.
