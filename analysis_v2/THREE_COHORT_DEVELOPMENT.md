# Mixed three-cohort development analysis

Complete strict-production Yachida profiles may be combined with historical
Feng and Zeller profiles to exercise the full downstream analysis before the
two new CRC-cohort production runs finish. This mixed-provenance engineering
run is always marked `DEVELOPMENT_ONLY`; it is not manuscript evidence.

The preparation runner rebuilds Yachida canonical evidence from its sealed
profiles and exact spike-design files, combines it with the validated legacy
Feng/Zeller canonical input, validates the union, and derives paired endpoints.
The existing map-reduce workflow then refits all models and reports.

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration

export LEGACY_INPUT_ROOT="$PWD/work/analysis_v2_legacy_native_complete_dev_20260909_171148"
export YACHIDA_ENV="$PWD/config/yachida.strict-production.env"
export YACHIDA_MANIFEST="$PWD/work/yachida_67x3/metadata/pilot_batched.tsv"
export YACHIDA_INDEPENDENT_MANIFEST="$PWD/work/yachida_67x3/metadata/independent_10_per_condition.tsv"
export THREE_COHORT_INPUT="$PWD/work/analysis_v2_three_cohort_input_dev_$(date +%Y%m%d_%H%M%S)"
export ANALYSIS_SIF=/mnt/beegfs/apptainer/images/ground_truth_analysis_v1.sif

OUTDIR="$THREE_COHORT_INPUT" \
  bash analysis_v2/prepare_three_cohort_development_input.sh

export DEV_INPUT_ROOT="$THREE_COHORT_INPUT/combined"
export OUTDIR="$PWD/work/analysis_v2_three_cohort_mapreduce_dev_$(date +%Y%m%d_%H%M%S)"

bash analysis_v2/submit_legacy_biomarker_mapreduce.sh
```

Keep `YACHIDA_MANIFEST` exported during submission so setup can build covariate
metadata for all three cohorts. The final seal depends on every report.
