# Parallel DA engineering pilot

Not definitive results or a statistical false-positive calibration. Prepare first,
inspect, then explicitly submit the array. No automatic production launch.

## Scope

84 independent model contexts, one CPU per fit, existing pinned MaAsLin2 wrapper:

- 12 DA1: three cohorts, two profilers, Adenoma/Control and CRC/Control, age/sex.
- 12 DA2: each cohort/profiler, adenoma background, full community dose .001 and
  ten individual Fnuc pairs at .0001. Both are nominal target dose .0001, not
  necessarily exactly equal achieved doses. Random sample intercept.
- 60 allocation contexts: 20 selected cohort/profiler/background allocations,
  each with baseline-only null, community total .0001 and .005, n=5/group.
  Strata cycle over cohort/background/profiler; first draw in all twelve strata,
  second in the first eight. Selection depends on no fitted response. Assignments
  come unchanged from the frozen full-pool allocation ledger; shared draws across
  profilers. The twenty contexts are not twenty independent biological cohorts.

Families are comparison-specific DA1 and full-background DA2/DA3. Abundance is
native fraction, missing species in valid sparse profiles are zero. No TSS or
recovery calibration. Full family BH and non-estimability policies are unchanged.
Warnings, singular/nonconverged paired fits or missing results block the context;
do not silently fall back to a different model. Review blocked diagnostics before
extending the pilot. With thousands of sparse species, paired failures are a
real possibility this pilot is intended to reveal.

## Preparation and dose records

prepare_da_pilot.sbatch consumes verified inventory and DA1 checkpoint roots,
rechecks all recorded source hashes and writes contexts, observations, families,
achieved_doses and nominal_matches, plus source hashes/checksums. Dose records
are canonical design-derived fractions, not independently measured read-origin
counts. Nominal community target matching divides total dose by ten (equal-panel
design); achieved values are preserved separately with differences. It is not an
exact achieved-dose matched comparison. Sensitivity matching policy and scientific
individual/community comparison analysis remain future work.

Environment: PROJECT, ANALYSIS_SIF, DA_INVENTORY_ROOT, DA1_CHECKPOINT_ROOT,
fresh DA_PILOT_PLAN. Preparation records source commit and container SHA.

## Parallel execution and resumption

After plan inspection: sbatch --array=0-83%4 --export=ALL analysis_v2/run_da_pilot.sbatch.
Export DA_PILOT_RESULTS as a fresh shared output root. Concurrency is configurable.
Threads are capped at one, no nested fitting parallelism. Default 16 GB/4 hours
per task are safety requests; measure actual resources with sacct, do not assume
they are final production needs. No short-fit batching yet: profile runtime first.

Each task verifies plan checksums, source profile hashes, image SHA, pinned package
version, and writes into a unique attempts directory. A per-context file lock
prevents duplicate fits. On success checksums cover model/input/diagnostic files,
then the complete directory is atomically published. Re-running the same index
verifies a completed result and checks matching code/plan provenance before reuse.
Failures retain FAILED.json, model logs and any backend diagnostics; no result
SUCCESS. Input/code changes require fresh outputs, not overwriting old results.

No final study seal or completeness/false-positive summary is generated yet.
Aggregation, complete family accounting, richer model diagnostics and review of
null rejection rates are the next steps after the pilot. Twenty selected null
contexts cannot establish error control; larger null validation still required.

## Validation

test_da_pilot.py uses a synthetic single-cohort fixture and mocked model runner;
checks disjoint groups, baseline-only nulls, dose-match records, resume and retained
failure diagnostics. It does not claim actual backend validation. Existing
test_maaslin_context.R runs the actual pinned backend when REQUIRE_MAASLIN_BACKEND=1.
Cluster container tests and real pilot results remain required. Model manifest
and array launcher have not yet been exercised against cluster study data.
