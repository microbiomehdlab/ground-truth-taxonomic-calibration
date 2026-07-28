# Public CRC cohorts

The spike-in benchmark uses two public paired-end stool-metagenome cohorts:

| Directory | Study | Samples used |
|---|---|---:|
| `fengq/` | Feng et al. 2015 | 154 |
| `zellerg/` | Zeller et al. 2014 | 156 |

Together these provide 310 biological backgrounds for the community-spike
experiment. The independent experiment uses a balanced subset of 60 samples
(30 per cohort), stratified across Control, Adenoma, and CRC backgrounds.

Each cohort directory contains:

- the study metadata used to associate sample IDs with diagnostic background;
- a downloader that organizes files by diagnostic group; and
- the generated EBI download-command inventory.

Review the original study access conditions and repository terms before
downloading or redistributing data.

## Download

```bash
bash datasets/fengq/fengq_download.sh /path/to/raw
bash datasets/zellerg/zellerg_download.sh /path/to/raw
```

Existing files are skipped. The resulting layout is:

```text
/path/to/raw/
  fengq/<condition>/<sample>/*_1.fastq.gz
  fengq/<condition>/<sample>/*_2.fastq.gz
  zellerg/<condition>/<sample>/*_1.fastq.gz
  zellerg/<condition>/<sample>/*_2.fastq.gz
```

## Build a profiler manifest

The dataset-agnostic manifest builder detects common paired-end naming
conventions:

```bash
bash spikes/scripts/spikein/make_samples_tsv.sh \
  /path/to/raw/fengq \
  /path/to/raw/zellerg \
  --subdir . \
  --maxdepth 4 \
  --id-depth 1 \
  > samples.tsv
```

Inspect the manifest and confirm every row has the expected pair before
submitting jobs:

```bash
head samples.tsv
wc -l samples.tsv
while IFS=$'\t' read -r sample r1 r2 rest; do
  [[ "$sample" == "sample_id" ]] && continue
  [[ -s "$r1" && -s "$r2" ]] || echo "Missing pair: $sample"
done < samples.tsv
```

For a publication run, archive the exact sample manifest, study metadata,
download date, source URLs, and checksums for all downloaded FASTQs.
