# Feng and Zeller definitive-analysis runbook

The Feng and Zeller implementations share one canonical-input builder,
upstream seal contract, readiness gate, and definitive downstream driver. The
cohorts remain separate analytical studies; shared code does not pool their
samples.

The frozen production sizes are 154 Feng samples and 156 Zeller samples. Each
cohort has a nested, outcome-independent 30-sample individual-spike subset
balanced as 10 Control, 10 Adenoma, and 10 CRC. Community spikes use every
eligible production sample.

The public upstream workflow for all three cohorts is
`UNIFIED_UPSTREAM_SEAL.md`, run through
`analysis_v2/run_cohort_upstream_audit.sbatch`. It consumes the native seal
created below as immutable provenance and emits the canonical
`production_seal_v2` layout that the evidence package and definitive downstream
workflows take as input. The native audit below is still how that source seal
is produced, and it remains unchanged. The frozen condition counts the unified
workflow requires are `Control=61,Adenoma=47,CRC=46` for Feng's 154 samples and
`Control=61,Adenoma=42,CRC=53` for Zeller's 156.

After a cohort finishes upstream, create its seal in a compute job because the
audit rehashes every retained output. For Feng:

```bash
export PROJECT="$PWD"
export CRC_ENV="$PWD/config/feng.strict-production.env"
export COHORT=feng
export PRODUCTION_MANIFEST="$PWD/datasets/fengq/manifests/production_manifest.tsv"
export INDEPENDENT_MANIFEST="$PWD/datasets/fengq/manifests/production_manifest.independent.tsv"
export EXPECTED_SAMPLES=154
AUDIT_JOB="$(sbatch --parsable --export=ALL run_crc_cohort_upstream_audit.sbatch)"
echo "$AUDIT_JOB"
```

For Zeller, set `COHORT=zeller`, use the `datasets/zellerg` manifests,
`EXPECTED_SAMPLES=156`, and the Zeller environment. The seal fails if any
sample lacks verified provenance, its retained-output receipt, top-level
success, or the exact native-profile topology: one baseline and seven community
profiles for every sample, plus 60 independent profiles for each member of the
frozen 30-sample subset. The audit also requires each retained
`sample_completion.tsv` to agree with the manifest identity, condition, subset
membership, design-row counts, and expected profile total. Production `SUCCESS`
sentinels may be zero-byte files created with `touch`; evidence tables must be
nonempty. The seal independently rehashes every file in every retained-output
receipt, verifies its byte count, and rejects any retained path inside the
sample's disposable scratch directory or outside that sample's persistent
results and QC roots. It invalidates an older seal before starting and writes
the complete sample-flow and covariate-missingness ledgers even when
incomplete.

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
