# Same-grid reference extension: all 100 saved allocations

Authorized scope: same three cohorts, both profilers, adenoma backgrounds,
uniform exposure U, n=10/20 per artificial group, nominal total community doses
0.0001/0.001/0.01, original conversion multiplier 1.0. Use all 100 original saved
allocations, with no new sample selection, random seed, filtering or DA method.
The completed pilot and half/double sensitivity remain separate evidence.

3600 reference contexts, 36000 target comparisons, 360 ten-context array tasks.
Preparation and collector bring the task total to 362. The original 20 allocations
are included in the full fresh run; their reference fits are recomputed under the
same method rather than relabeling outputs with another plan identity. Observed
primary fits are reused in every context. No FASTQs or taxonomic profiles rerun.

The shared launcher now accepts REFERENCE_ALLOCATIONS=20 (backward-compatible
default) or 100. Exposure remains fixed U. It derives array length from the fixed
declared grid; preparation independently checks generated tasks/context count
before releasing the array. Every worker snapshots the calculation code and
retains image/plan hash checks and atomic resume. Existing output roots are refused.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration
git pull --ff-only origin revised-analysis-v2
export PROJECT="$PWD"
export ANALYSIS_SIF="/mnt/beegfs/apptainer/images/ground_truth_analysis_v2_1.sif"
export DA3_BATCH_PLAN="$PWD/work/preproduction_checks_20261002T192903Z/draft_da3_plan"
export DA_INVENTORY_ROOT="$PWD/work/three_cohort_da_inventory_20261001T151231Z"
export GEFF_ROOT="$PWD/work/three_cohort_sealed_geff_20261001T135712Z"
export DA3_OBSERVED_REPORT="$PWD/work/da3_direct_maaslin_20261003T170501Z/REPORT_20261003T170503Z"
export REFERENCE_ALLOCATIONS=100
export REFERENCE_CONCURRENCY=72
export REFERENCE_ROOT="$PWD/work/reference_response_100_$(date -u +%Y%m%dT%H%M%SZ)"
bash analysis_v2/submit_reference_extension.sh
```

Concurrency is tasks, not dedicated nodes; each worker requests 1 CPU/12 GB.
Use only an allowed cluster budget. Dependencies require successful preparation
and every array task. Do not remove dependencies or overwrite failed outputs.

After the collector finishes:

```bash
cat "$REFERENCE_ROOT.jobs.tsv"
cat "$REFERENCE_ROOT/REPORT/status.json"
(cd "$REFERENCE_ROOT/REPORT" && sha256sum -c --quiet SHA256SUMS)
printf '%s/REPORT\n' "$REFERENCE_ROOT"
```

Expected status has contexts=3600, target_rows=36000. Download only REPORT for
local scientific review. This is conditional sample-selection stability, not
100 independent biological replications or population clinical power. Retain all
taxa, cohort differences, coefficients/SE/raw p/q, HC3 and separate NE states.
No automatic claim update, heterogeneous exposure expansion or production flag.

Local verification: complete 20- and 100-allocation synthetic prepare/run/resume/
collect pipelines (R process mocked, not actual MaAsLin), exact selected-grid
coverage and preservation of pilot allocations, mocked Slurm submissions for
72/360 tasks and dependency ordering, and invalid-allocation refusal before any
submission. Actual pinned backend remains a required cluster preparation fixture.
