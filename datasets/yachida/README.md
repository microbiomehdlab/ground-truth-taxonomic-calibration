# Yachida 2019 external reproducibility cohort

This workflow uses the 80 paired-end stool metagenomes deposited under
`DRA006684` for Yachida et al. (2019): 40 healthy controls and 40 stage III/IV
CRC cases. This deposited WGS subset has no adenoma group. It is therefore a
two-background external validation, not an exact repetition of the original
Control/Adenoma/CRC design.

## 1. Freeze the complete cohort before downloading reads

```bash
mkdir -p work/yachida/metadata
python3 datasets/yachida/build_manifest.py \
  --cache-dir work/yachida/metadata/source \
  --output work/yachida/metadata/cohort_manifest.tsv
```

The builder retrieves the DDBJ run, experiment, and sample XML plus the ENA
FASTQ inventory. It requires exactly 80 runs and the expected 40/40 condition
split, selects the paired `_1`/`_2` FASTQs, and writes a manifest SHA-256 file.
Archive this manifest with the run. Never rebuild sample selection from only
the files that happen to fit on disk.

## 2. Freeze the independent-spike subset

```bash
python3 scripts/select_samples_deterministically.py \
  --manifest work/yachida/metadata/cohort_manifest.tsv \
  --output work/yachida/metadata/independent_10_per_condition.tsv \
  --per-condition 10 \
  --selection-seed ground-truth-taxonomic-calibration-yachida-v1
```

Selection is ranked by a SHA-256 hash of the frozen seed and run accession. It
is independent of input row order, batch size, download order, and currently
available files. The community design uses all 80 manifest rows.

## 3. Use a verified streaming lifecycle

Process one sample (or a deliberately small batch) completely:

1. download both FASTQs and verify archive-provided MD5 and byte counts;
2. run MetaPrep;
3. profile the unspiked sample with both profilers;
4. if selected for the independent design, generate/profile one target at a
   time and remove its verified spiked FASTQs before generating the next;
5. generate/profile the community fractions;
6. retain compact profiler tables, logs, tool/database versions, spike-design
   files, and checksums outside disposable scratch;
7. delete raw/cleaned/spiked FASTQs only after all receipt entries verify.

`stream_sample.py` enforces the download and cleanup boundary:

```bash
python3 datasets/yachida/stream_sample.py \
  --manifest work/yachida/metadata/cohort_manifest.tsv \
  --sample-id DRR127476 \
  --scratch-root /path/to/disposable/yachida_scratch \
  --state-dir /path/to/persistent/yachida_state \
  --runner /path/to/site_yachida_runner.sh \
  --delete-inputs-after-verification
```

The runner is site-specific because MetaPrep is maintained in a separate
repository. It receives these environment variables:

```text
SAMPLE_ID RAW_R1 RAW_R2 SAMPLE_WORK RECEIPT TARGET_CONDITION STUDY
```

It must block until preprocessing, spike generation, both profiler runs, and
post-processing have completed. It then writes `RECEIPT` with:

```text
path    sha256    bytes
/persistent/path/to/output.tsv    <sha256>    <byte-count>
```

Every listed path must be non-empty and outside `SAMPLE_WORK`. The streaming
driver independently verifies all hashes and sizes. Cleanup is opt-in and is
refused without the per-sample safety sentinel. Omit
`--delete-inputs-after-verification` for the first test sample.

## Seed invariance

Spike-pool ART seeds and seqtk subsampling seeds use `stable-seed-v1`. They are
derived from `SEED_BASE` and stable identifiers (sample, target/community,
fraction, and community member). Slurm task IDs, manifest rows, batch numbers,
fraction-list positions, and local paths are deliberately excluded. Thus the
same sample produces the same spike reads whether it is processed alone,
first, last, or after a resumed run.

## Recommended validation order

Start with one control and one CRC sample without automatic cleanup. Compare a
second run made from reordered two-row manifests. Spike-design seeds and output
checksums must match. Then enable verified cleanup and expand to the frozen
20-sample independent subset and 80-sample community cohort.
