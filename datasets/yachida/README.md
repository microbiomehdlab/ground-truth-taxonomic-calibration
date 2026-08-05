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
more expensive independent-species design.

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
`43a80d6e22f1a1ffd33e9d8e6c58ab9d37ad7b0c`. It records configurable output
placement and deterministic discovery of paired host-depleted reads. The
adapter requires this exact revision and refuses uncommitted changes to tracked
MetaShotgunPrep files; untracked interpreter caches do not affect validation.

`run_metashotgunprep.sh` adapts the streaming R1/R2 inputs to MetaShotgunPrep's
required directory layout without duplicating the downloads. It verifies the
pinned MetaShotgunPrep revision and every Bowtie2 index component, deliberately
does not request MetaShotgunPrep raw deletion, validates both cleaned gzip
mates, preserves compact FastQC/log provenance, and writes
`metashotgunprep_outputs.env` for downstream profiling.

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

## Seed invariance

ART and seqtk use `stable-seed-v1`, derived from stable sample, target,
community, fraction, and member identifiers. Slurm task IDs, manifest rows,
batch numbers, fraction positions, and local paths are excluded. The same
sample therefore produces the same spike reads when run alone, in another
batch, or after resumption.
