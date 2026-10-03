# Full direct-MaAsLin2 DA3 run

Approved by the user on 3 October 2026 after the direct-package timing trial:
48/48 contexts complete, zero mismatched features. That trial projected about
774 single-worker hours / 26 ideal hours at 30 workers and 724 GiB with full
model/input retention. Those are estimates, not scheduler guarantees.

This runner executes the **same 120,000 DA3 experiments** through actual
MaAsLin2 1.18.0. It does not regenerate allocations, change features/doses,
rerun profiling, replace the old exploratory results, or repeat permutations.
These are artificial independent-person comparisons in adenoma/CRC backgrounds,
not actual disease-association rates. Clinical DA1 and paired DA2 are unchanged.

## Execution and retention

- Copy the existing checksummed 4,800-batch plan, preserving every context.
- Snapshot committed source and pin/check the analysis image and source hashes.
- One R process per batch of 25 contexts amortizes package startup. Six Slurm
  arrays of 800 tasks share the requested total concurrency (default 30).
- One CPU / 8 GiB / two-hour limit per batch; BLAS/OpenMP threads are one.
- The exact timing-pilot wrapper is reused: explicit `log2(1+a/1e-8)`, group-only
  unpaired LM, no additional normalization/filtering, full-frozen-family BH.
  Native package q-values remain separately labeled. Constant/perfect fixed-fit
  raw p is NA; multiplicity bookkeeping p is 1. These are not package defaults.
- Models are temporarily saved and inspected **exactly as in the timing pilot**.
  Successful-context native models and duplicate native tables are then removed
  from node-local scratch. This changes retention, not fitting or diagnostics.
- Persist gzip-compressed statistics for **all family members**, target tables,
  comparisons that disagree with the fast LM, per-context status/count/timing,
  compressed R console logs, session versions, input hashes and run identity.
  Failure diagnostics and inputs are preserved in `attempts/`; no failed batch
  gets promoted to `results/`.
- Successful outputs are atomically promoted and locked per batch. Resume checks
  hashes, identity, coverage and source inputs before skipping completed work.
  A corrupt completed batch is refused, never silently overwritten.

Model/storage savings have not yet been measured on the full run; do not reuse
the old 724 GiB estimate as the compact output prediction. Temporary models
still consume node scratch, and failure evidence may be larger.

## Submit once (André runs this on the cluster)

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA3_BATCH_PLAN="$PWD/work/preproduction_checks_20261002T192903Z/draft_da3_plan"
export MAASLIN_RUN_ROOT="$PWD/work/da3_direct_maaslin_$(date -u +%Y%m%dT%H%M%SZ)"
export MAASLIN_CONCURRENCY=30
bash analysis_v2/submit_maaslin_production.sh
```

The preparation job validates the copied plan and runs a mandatory actual-package
fixture confirming compact batched results equal the original wrapper. Arrays
depend on successful preparation. This is an operational gate, not another
statistical-method campaign. If it fails, a `PREP_FAILURE/README.txt` points to
its log; arrays remain dependency-held and no large fit runs. Do not manually
release dependencies. Report collection depends on **any** array outcome so a
failed batch is explicitly reported rather than counted as a nondiscovery.

`jobs.tsv` records every submitted ID. `report_paths.txt` records versioned report
directories; shell variables are not persistent across sessions. After jobs finish:

```bash
REPORT_ROOT="$(tail -n 1 "$MAASLIN_RUN_ROOT/report_paths.txt")"
cat "$REPORT_ROOT/status.json"
(cd "$REPORT_ROOT" && sha256sum -c --quiet SHA256SUMS)
```

A complete run must report 4,800 batches and 120,000 contexts, no failures;
any numerical/significance disagreement is flagged separately for review.
`production_authorized: false` in this final scientific-review report does not
mean the user-approved computation was withheld: it means completed computation
does not automatically license manuscript conclusions or a new analysis.

## Resume

After all jobs for this root have stopped, restore its exact absolute path and
the original image, then use:

```bash
bash "$MAASLIN_RUN_ROOT/source/analysis_v2/submit_maaslin_production.sh" --resume
```

Resume uses the **original source snapshot and copied plan**, not a new checkout.
It queues the batch grid, skips verified completed batches without fitting again,
and retries failed/unstarted batches in new attempts. Failed attempts and previous
reports are retained. It refuses submission if recorded jobs are still active.
Preparation failure before a valid identity exists requires a new root after
diagnosing the failure; do not alter an old identity to force continuation.

## Outputs and review

Each `results/batch_XXXXX/features.tsv.gz` retains feature names, coefficient,
standard error, raw p, native q, estimability/status, BH bookkeeping p,
full-family q and positive-discovery indicator with exact context ID.
The plan supplies all metadata and sample allocations. `targets.tsv.gz` adds
the ten target labels. `differences.tsv.gz` contains only implementation
disagreements; it has a header even when empty. Logs and session metadata
preserve diagnostic evidence without retaining all model objects.

The consolidated report provides `context_summary.tsv`, `targets.tsv.gz`,
`feature_results_inventory.tsv` (checksummed full-feature batch files), identity,
status and checksums. Full-family tables remain sharded to avoid a duplicate
giant table. Partial reports are labeled INCOMPLETE and must not be used for
rates. Keep the whole run root; downloading only REPORT is sufficient for an
initial review, not for archiving complete feature-level results.

Local R fixtures test numerical equivalence through a test backend and retention;
Python fixtures test failure preservation, immutable inputs, corrupted-output
rejection, resume and collection. Actual package execution is mandatory on the
cluster. Existing small-n/heteroscedasticity limitations remain reportable;
direct MaAsLin2 use does not resolve those statistical limitations.
