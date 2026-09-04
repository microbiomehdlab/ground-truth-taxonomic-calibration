# Yachida assembly-sensitivity runbook

This workflow consumes immutable upstream evidence and writes a new downstream
directory. It does not alter the original or clean-arm profiles. The semantic
audit that passed 720 clean-arm native profiles is an important development
gate; the definitive run repeats that audit on the exact profiles represented
in its canonical table.

On Lobo, after the strict original Yachida experiment is complete and sealed:

```bash
cd /mnt/nfs/microbiomehd/crc-lab/projects/ground-truth-taxonomic-calibration

export SENSITIVITY_ROOT="$PWD/work/yachida_assembly_sensitivity_20260901"
export ORIGINAL_ROOT="$PWD/work/yachida_strict_final_20260823"
export BASELINE_ROOT="$ORIGINAL_ROOT"
export ANALYSIS_SIF=/mnt/beegfs/apptainer/images/ground_truth_analysis_v1.sif
export OUTDIR="$PWD/work/analysis_v2_yachida_assembly_sensitivity_$(date +%Y%m%d_%H%M%S)"

bash analysis_v2/run_assembly_sensitivity.sh
```

Do not reuse `work/analysis_v2_development_20260904` as the output directory.
That directory is explicitly development-only. A successful definitive run
contains top-level `SUCCESS`, canonical validation, an exact native-profile
semantics audit, paired endpoints, model diagnostics and estimates, provenance,
and SHA-256 seals. Failure at any gate prevents the final `SUCCESS` marker.

Before interpreting estimates, verify:

```bash
test -s "$OUTDIR/SUCCESS"
test -s "$OUTDIR/canonical/validation/SUCCESS"
test -s "$OUTDIR/endpoints/SUCCESS"
test -s "$OUTDIR/models/assembly_sensitivity/SUCCESS"
grep -P '^status\tPASS$' "$OUTDIR/provenance/run_manifest.tsv"
```

The primary table is
`models/assembly_sensitivity/assembly_slope_contrasts.tsv`. The two
`clean_minus_original` rows quantify assembly-choice sensitivity separately for
each profiler. They do not isolate contamination from strain or database
representation effects.
