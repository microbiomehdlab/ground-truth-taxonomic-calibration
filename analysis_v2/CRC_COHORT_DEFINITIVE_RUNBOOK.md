# Feng and Zeller definitive-analysis runbook

The Feng and Zeller implementations share one canonical-input builder,
upstream seal contract, readiness gate, and definitive downstream driver. The
cohorts remain separate analytical studies; shared code does not pool their
samples.

The frozen production sizes are 154 Feng samples and 156 Zeller samples. Each
cohort has a nested, outcome-independent 30-sample individual-spike subset
balanced as 10 Control, 10 Adenoma, and 10 CRC. Community spikes use every
eligible production sample.

After a cohort finishes upstream, create its seal. For Feng:

```bash
source config/feng.strict-production.env
python3 analysis_v2/scripts/seal_crc_cohort_upstream.py \
  --cohort feng \
  --manifest datasets/fengq/manifests/production_manifest.tsv \
  --independent-manifest datasets/fengq/manifests/production_manifest.independent.tsv \
  --state-dir "$CRC_STATE_DIR" \
  --results-root "$PERSISTENT_RESULTS_ROOT" \
  --expected-samples 154 \
  --outdir "$CRC_STATE_DIR/production_seal"
```

For Zeller, use `--cohort zeller`, the `datasets/zellerg` manifests,
`--expected-samples 156`, and the Zeller environment. The seal fails if any
sample lacks verified provenance, its retained-output receipt, top-level
success, or exactly 8/68 successful native profiles. It writes the complete
sample-flow and covariate-missingness ledgers even when incomplete.

Then run a preflight in a new directory:

```bash
export PROJECT="$PWD"
export COHORT=feng
export CRC_ENV="$PWD/config/feng.strict-production.env"
export ANALYSIS_SIF=/mnt/beegfs/apptainer/images/ground_truth_analysis_v1.sif
export RUN_ROOT="$PWD/work/analysis_v2_feng_preflight_$(date +%Y%m%d_%H%M%S)"
export PREFLIGHT_ONLY=1
sbatch --export=ALL run_crc_cohort_definitive_analysis.sbatch
```

Use the corresponding Zeller values for that cohort. After preflight passes,
unset `PREFLIGHT_ONLY`, select a fresh `RUN_ROOT`, and submit again. Never turn
a partial or development run into a definitive package.

The definitive driver builds and validates the canonical table, audits exact
native profiles, derives paired endpoints, runs independent and community
detection/continuous models, performs pooled-primary and phenotype-stratified
artificial-biomarker analyses, evaluates disease-biomarker propagation,
generates reports and calibration linkage, and creates a checksum seal.

Once all three cohort packages exist, set `YACHIDA_RUN`, `FENG_RUN`,
`ZELLER_RUN`, `ANALYSIS_SIF`, and a fresh `OUTDIR`, then execute
`analysis_v2/run_three_cohort_definitive_synthesis.sh`. It refuses development
or incomplete cohort packages.
