# Yachida 2019 external reproducibility cohort

The complete sequencing project `PRJDB4176` contains 645 runs. The curated
Segata metadata file `Public_study__YachidaS_2019.tsv` defines the 615 samples
used here: 290 Control, 67 Adenoma, and 258 CRC. The remaining 30 project runs
are not silently included; the manifest builder records them separately.

The initial external-validation pilot is balanced at the size of the smallest
group: 67 Control, all 67 Adenoma, and 67 CRC samples (201 total).

## 1. Freeze the complete curated cohort

Copy `Public_study__YachidaS_2019.tsv` to the cluster, then run:

```bash
mkdir -p work/yachida/metadata

python3 datasets/yachida/build_manifest.py \
  --study-metadata /path/to/Public_study__YachidaS_2019.tsv \
  --cache-dir work/yachida/metadata/source \
  --output work/yachida/metadata/cohort_manifest.tsv
```

The builder downloads the official ENA `PRJDB4176` FASTQ inventory, requires a
one-to-one run mapping for every curated BioSample, selects paired `_1`/`_2`
FASTQs, validates the 290/67/258 counts, and writes:

- `cohort_manifest.tsv` and its SHA-256 file;
- `cohort_manifest.excluded_project_runs.tsv`, containing the 30 project runs
  absent from the curated metadata, and its SHA-256 file;
- `cohort_manifest.provenance.tsv`, recording the source metadata and ENA
  inventory checksums, schema version, and verified counts.

Archive these files with the run. Never select samples from only the FASTQs
that happen to be present on disk.

## 2. Freeze the matched 67/67/67 pilot

```bash
python3 datasets/yachida/build_pilot_design.py \
  --manifest work/yachida/metadata/cohort_manifest.tsv \
  --output work/yachida/metadata/pilot_67_per_condition.tsv \
  --per-condition 67 \
  --reference-condition Adenoma \
  --selection-seed ground-truth-taxonomic-calibration-yachida-pilot-v1
```

All adenomas are retained. Control and CRC samples are selected so their exact
sex × age-bin counts match the adenoma group. Fixed age bins are `<50`,
`50–59`, `60–69`, and `≥70`. Within each stratum, a stable accession hash
determines selection. Selection is therefore invariant to metadata row order,
download order, batch size, and available disk contents.

## 3. Freeze the nested independent-spike subset

```bash
python3 scripts/select_samples_deterministically.py \
  --manifest work/yachida/metadata/pilot_67_per_condition.tsv \
  --output work/yachida/metadata/independent_10_per_condition.tsv \
  --per-condition 10 \
  --selection-seed ground-truth-taxonomic-calibration-yachida-independent-v1
```

This selects 10 Control, 10 Adenoma, and 10 CRC samples from within the frozen
201-sample pilot. The recommended first analysis uses all 201 samples for
baseline and community-spike profiling, and these nested 30 samples for the
more expensive independent-species design. This second selection is balanced
by condition but is not rematched on age or sex. It uses no profiling,
recovery, or biomarker outcome.

Before production, freeze both selection stages in an unambiguous provenance
table and inspect descriptive covariate balance:

```bash
python3 datasets/yachida/audit_selection_balance.py \
  --pilot-manifest work/yachida_67x3/metadata/pilot_67_per_condition.tsv \
  --independent-manifest work/yachida_67x3/metadata/independent_10_per_condition.tsv \
  --output-dir work/yachida_67x3/selection_audit

cd work/yachida_67x3/selection_audit
sha256sum -c SHA256SUMS
```

The generated report is descriptive rather than an outcome-dependent pass/fail
test. Decide whether the frozen subset is acceptable before examining recovery
or biomarker results. The output is private run evidence under Git-ignored
`work/`.

## 4. Create processing batches of at most 10

Use a new, empty batch directory:

```bash
python3 scripts/assign_processing_batches.py \
  --manifest work/yachida/metadata/pilot_67_per_condition.tsv \
  --output work/yachida/metadata/pilot_batched.tsv \
  --batch-dir work/yachida/metadata/batches \
  --max-batch-size 10 \
  --batch-seed ground-truth-taxonomic-calibration-yachida-batches-v1
```

For 201 samples this produces 21 deterministically assigned batches: 12
batches of 10 and 9 batches of 9. Conditions are round-robin interleaved. Batch
assignment never enters sample selection or spike-seed derivation, so failed
samples may be rerun individually without changing their synthetic reads.

## 5. Use a verified streaming lifecycle

Create the local, ignored site configuration:

```bash
cp datasets/yachida/yachida.env.example config/yachida.env
```

Edit `config/yachida.env` with the local MetaShotgunPrep checkout, pinned Git
revision, Bowtie2 host-index prefix, and persistent QC directory. Repository
documentation and committed scripts use placeholders only; host-specific paths
remain exclusively in this ignored file.

The pinned MetaShotgunPrep revision is
`5923619824799457e89bd78b211ed481b7cb6f3f`. It includes configurable output
placement, deterministic discovery of paired host-depleted reads, and
standardized persistent QC exports. The
adapter requires this exact revision and refuses uncommitted changes to tracked
MetaShotgunPrep files; untracked interpreter caches do not affect validation.

`run_metashotgunprep.sh` adapts the streaming R1/R2 inputs to MetaShotgunPrep's
required directory layout without duplicating the downloads. It verifies the
pinned MetaShotgunPrep revision and every Bowtie2 index component, deliberately
does not request MetaShotgunPrep raw deletion, performs a complete
FASTQ-structure and normalized mate-identifier synchronization scan, preserves
the resulting pair count/read-characteristics table with compact FastQC/log
provenance, and writes
`metashotgunprep_outputs.env` for downstream profiling.

The default technical eligibility rule is intentionally non-selective:
`MIN_POST_QC_PAIRS=1` retains every non-empty, structurally valid synchronized
library. If a higher sequencing-depth exclusion is scientifically required,
set and document it before processing any cohort; do not choose it after
examining biomarker results.

The frozen Yachida profiling protocol uses `BRACKEN_THRESHOLD=10`. This is
Bracken's minimum read-count threshold (`-t`), not a thread count. Kraken2 CPU
allocation is controlled independently by `K2_THREADS`. The pre-production
gate and per-profile runner both fail if the threshold differs from 10, and
every profile records it in `profiling_parameters.tsv`.

For every sample in the current batch:

1. download both FASTQs and verify archive-provided MD5 and byte counts;
2. run metaprep;
3. profile the unspiked sample with both profilers;
4. for the nested independent subset, generate/profile one target at a time
   and remove its verified spiked FASTQs before generating the next target;
5. generate/profile the community fractions;
6. retain compact profiler tables, logs, tool/database versions, spike-design
   files, and checksums outside disposable scratch;
7. delete raw, cleaned, and spiked FASTQs only after all receipt entries verify.

`stream_sample.py` enforces the download and cleanup boundary:

```bash
python3 datasets/yachida/stream_sample.py \
  --manifest work/yachida/metadata/pilot_batched.tsv \
  --sample-id SAMD00115014 \
  --scratch-root /path/to/disposable/yachida_scratch \
  --state-dir /path/to/persistent/yachida_state \
  --runner /path/to/site_yachida_runner.sh \
  --delete-inputs-after-verification
```

The site-specific runner is required because metaprep is maintained in a
separate repository. It receives:

```text
SAMPLE_ID RAW_R1 RAW_R2 SAMPLE_WORK RECEIPT TARGET_CONDITION STUDY
```

It must block until processing finishes and write `RECEIPT` with header
`path<TAB>sha256<TAB>bytes`, listing every persistent output. The driver
independently verifies these files. Cleanup is opt-in, restricted to the exact
per-sample scratch directory, and refused without its safety sentinel. Omit
`--delete-inputs-after-verification` for the first validation batch.

### Preprocessing smoke test

Before merging a new MetaShotgunPrep revision or running profiling, submit one
selected sample through the non-deleting smoke-test runner:

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
export SAMPLE_ID=SAMD00164833
sbatch --export=ALL run_yachida_preprocessing_smoke.sbatch
```

The job uses the source-built upstream image, verifies the manifest download,
checks both cleaned gzip mates and MetaShotgunPrep's completion marker, retains
QC with checksums, and deliberately leaves raw and cleaned reads under
`work/yachida_67x3/smoke_scratch/` for inspection. It does not run profiling or
spike generation and cannot authorize deletion of its scratch directory.

After preprocessing passes, reuse the retained cleaned mates for a non-deleting
baseline profiler smoke test:

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
export SAMPLE_ID=SAMD00164833
sbatch --export=ALL run_yachida_baseline_smoke.sbatch
```

This runs Kraken2/Bracken and MetaPhlAn directly from the pinned upstream image,
uses the configured UHGG and vJan25 databases, validates species-level output,
and checksums the compact profiles under `work/yachida_67x3/baseline_smoke/`.
It does not regenerate preprocessing, create spikes, or delete cleaned reads.

Before the spike smoke test, freshly download and audit the exact versioned
target assemblies using the pinned NCBI Datasets CLI in the upstream image, as
documented in `spikes/README.md`. Spike pools are shared immutable assets and
must be generated once from those audited FASTAs using a new parameter-specific
pool directory.

Create the local Yachida spike configuration without editing a template:

```bash
bash spikes/scripts/spikein/configure_spike_workflow.sh \
  --site-env config/yachida.env \
  --work-dir work/yachida_67x3 \
  --reference-dir references/genomes
```

The generated `work/yachida_67x3/spikein.env` deliberately defers
`SAMPLES_TSV`: the cohort manifest contains download metadata, not cleaned-read
paths. The streaming runner supplies a three-column, batch-specific cleaned-read
manifest when it performs spike insertion. Pool generation does not require a
sample manifest.

After finalizing the shared pools, run the non-deleting end-to-end spike smoke
test on the retained preprocessing-smoke sample:

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
export SAMPLE_ID=SAMD00164833
sbatch --export=ALL run_yachida_spike_smoke.sbatch
```

The job generates an independent Fnuc spike and a ten-member community spike,
both at 0.01% total abundance, validates their paired FASTQs, and profiles each
with Kraken2/Bracken and MetaPhlAn 4. Synthetic FASTQs are retained for this
first audit. No cohort-scale processing or deletion is authorized by this job.
The production workflow uses deterministic paired single-pass selection and
records its sampling mode and seeds in the retained provenance.

Before cohort-scale use, run the supplied integration check on small temporary
background and pool subsets:

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
export SAMPLE_ID=SAMD00164833
sbatch --export=ALL run_yachida_sampling_integration.sbatch
```

This verifies deterministic paired selection, FASTQ construction, and design
records without modifying the finalized pools or completed smoke run.

### Required pre-production gate

After the smoke and integration tests, run the one-time fail-closed gate before
submitting any additional cohort batches:

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
export SPIKE_ENV="$PWD/work/yachida_67x3/spikein.env"
AUDIT_JOB=$(sbatch --parsable --export=ALL \
  run_yachida_preproduction_audit.sbatch)
sacct -j "$AUDIT_JOB" --format=JobID,JobName,State,Elapsed,ExitCode
```

The job verifies the frozen 67/67/67 pilot, nested 10/10/10 independent subset,
ten-target panel, synchronized pool-count index, checksum-manifest seal, all 20
pool file checksums, complete pool FASTQ/mate-ID integrity, and deterministic
exact community allocation for integer totals from 1 through 1000. It writes a
private, Git-ignored report and pool read-characteristics tables under
`work/yachida_67x3/preproduction_audit/`. A `COMPLETED` Slurm state, exit code
`0:0`, and `SUCCESS` file are required before production. Set
`VERIFY_POOL_CONTENTS=0` only for subsequent repeat audits after the same pool
files have already passed the full checksum scan.

## 6. Storage-aware full batch lifecycle

The full runner processes one sample per Slurm array task and one target at a
time. For each target, all currently missing fractions are selected in one
pool scan. A configurable number of fractions are then profiled concurrently;
each synthetic pair is deleted immediately after both profilers have produced
and independently verified their compact receipt. Community fractions follow
the same lifecycle. Raw and cleaned reads remain protected by
`stream_sample.py` until all expected sample outputs pass a final persistent
receipt.

First create the one-time pool pair-count index. This validates synchronized
mate counts in one scan and avoids rescanning the large finalized pools merely
to determine their sizes during every fraction. The finalization checksum
manifest must already exist; the index records its checksum without repeating
the separate full-file hash scan:

```bash
export PROJECT="$PWD"
export SPIKE_ENV="$PWD/work/yachida_67x3/spikein.env"
sbatch --export=ALL run_index_yachida_pools.sbatch
```

Configure these ignored site paths in `config/yachida.env`:

```text
PERSISTENT_RESULTS_ROOT=/private/persistent/yachida/results
YACHIDA_SCRATCH_ROOT=/disposable/yachida/scratch
YACHIDA_STATE_DIR=/private/persistent/yachida/state
```

Disposable scratch should be on the site's high-throughput scratch filesystem.
Persistent results, QC, and state should remain on backed-up storage. Finalized
pools may optionally be cached on the same fast filesystem when site capacity
permits. Stage and verify that immutable cache in a Slurm job:

```bash
export PROJECT="$PWD"
export SOURCE_POOLS_DIR=/path/to/finalized/spike_pools_readlen100_cov2000
export DESTINATION_POOLS_DIR=/path/to/fast/cache/spike_pools_readlen100_cov2000
sbatch --export=ALL run_stage_yachida_pools.sbatch
```

Use these production settings:

```text
YACHIDA_POOLS_DIR=/path/to/finalized/spike_pools_readlen100_cov2000
FASTQ_ASSEMBLY_MODE=gzip_members
PROFILE_CONCURRENCY=2
```

`gzip_members` produces a standards-compliant multi-member gzip file without a
redundant decompression/recompression pass. The two-way profile configuration
requests 32 CPUs and 96 GB per sample task. This allocation accommodates two
concurrent profiler processes while allowing tasks to run on cluster nodes with
approximately 126 GB RAM; it changes scheduling only, not analytical output.
The runner records scratch bytes
after preprocessing, construction, and verified cleanup in each sample's
`scratch_usage.tsv`; use the observed maximum plus at least 25% headroom when
choosing array concurrency. Begin with one task, then raise concurrency only
after confirming filesystem and scheduler capacity.

Validate the first position of a frozen batch without deleting raw or cleaned
reads:

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
bash datasets/yachida/submit_batch.sh \
  --manifest work/yachida_67x3/metadata/batches/batch_001.tsv \
  --array-indices 1 \
  --max-concurrent 1
```

After inspecting that sample, rerun the full batch with verified cleanup
enabled. The already verified first task is not recomputed; it only removes its
sentinel-protected scratch if it still exists:

```bash
bash datasets/yachida/submit_batch.sh \
  --manifest work/yachida_67x3/metadata/batches/batch_001.tsv \
  --max-concurrent 5 \
  --delete-verified-inputs
```

Only increase `--max-concurrent` after measuring peak scratch usage. Failed or
interrupted tasks retain their exact per-sample scratch and resume from verified
profile receipts. They never authorize deletion.

After every array task completes, independently seal the batch in a short
Slurm job (the audit rehashes compact retained outputs and should not be run on
a login node):

```bash
export PROJECT="$PWD"
export YACHIDA_ENV="$PWD/config/yachida.env"
export BATCH_MANIFEST="$PWD/work/yachida_67x3/metadata/batches/batch_001.tsv"
sbatch --export=ALL run_yachida_batch_audit.sbatch
```

The audit rehashes every retained output and writes a batch-level `SUCCESS`
receipt. A sample outside the nested 30-sample independent subset must have one
baseline and seven community profiles. A nested sample must additionally have
60 independent profiles (ten taxa by six fractions).

When a validation rule is added after early samples were completed, audit which
samples lack its persistent evidence rather than assuming uniform processing:

```bash
source config/yachida.env
python3 datasets/yachida/audit_processing_consistency.py \
  --manifest work/yachida_67x3/metadata/pilot_batched.tsv \
  --results-root "$PERSISTENT_RESULTS_ROOT" \
  --qc-root "$PERSISTENT_QC_ROOT" \
  --output work/yachida_67x3/processing_consistency.tsv
```

Rows marked `reprocess_with_final_pipeline` have completed results but lack the
final paired-FASTQ integrity record. Reprocess those samples from verified raw
inputs before forming the manuscript dataset. Do not mix them silently with
final-pipeline outputs.

### Optional target-abundance equivalence audit

To compare two independently generated result trees for the same nested sample,
use the read-only target audit. It resolves profiler-specific aliases, checks all
60 independent and 70 community target values per profiler, treats features
absent from a profiler table as zero abundance, and reports values in percent:

```bash
python3 datasets/yachida/compare_target_abundances.py \
  --reference-root /path/to/reference/results/YachidaS_2019/SAMPLE_ID \
  --candidate-root /path/to/candidate/results/YachidaS_2019/SAMPLE_ID \
  --panel spikes/spike_panel.tsv \
  --aliases examples/spike_taxon_aliases.csv \
  --output work/private_validation/SAMPLE_ID.target_abundance_comparison.tsv \
  --tolerance-pp 0.001
```

The generated comparison is local run output and is ignored by Git.

## Seed invariance

ART and seqtk use `stable-seed-v1`, derived from stable sample, target,
community, fraction, and member identifiers. Slurm task IDs, manifest rows,
batch numbers, fraction positions, and local paths are excluded. The same
sample therefore produces the same spike reads when run alone, in another
batch, or after resumption.
