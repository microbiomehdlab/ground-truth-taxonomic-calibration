# Three-cohort MetaPhlAn recovery-only derivation

Reference preparation passed on the cluster as job 3097829, 2026-10-01:
work/three_cohort_sealed_geff_20261001T135712Z, 511 samples, empty error log,
root checksums verified by André. Yachida 201, Feng 154, Zeller 156.

run_three_cohort_metaphlan_recovery.sbatch consumes this immutable bundle and
the definitive shared canonical preflights at timestamp 20260929T155304Z.
It checks source hashes and exact baseline identities against G_eff selection,
reference values and coverage, then derives endpoints under genome_equivalent.
No DA models, feature filtering, calibration or upstream profiling are run.

Both comparisons are present in each cohort's endpoints/paired_endpoints.tsv:

- MetaPhlAn primary: response_ratio_profiler_scale, with
  recovered_spike_signal_profiler_scale / implanted_signal_profiler_scale.
- Old sensitivity: response_ratio, using the unchanged read-fraction expectation.
- Neither comparison changes observed_native_abundance.

These tables also contain Bracken rows from the original read-proportional
endpoint calculation. They are NOT the new all-input Bracken primary results;
that measurement remains in the separately verified Bracken recovery bundle.
Select profiler == metaphlan4 for this report. Retain negative recovery ratios.

Required exports: PROJECT, ANALYSIS_SIF, GEFF_ROOT, PREFLIGHT_STAMP,
and a fresh METAPHLAN_RECOVERY_ROOT. All data paths must be visible under the
site wrapper's /mnt/nfs bind. Tests: test_geff_baseline_identity.py,
test_geff_canonical_integration.py and test_metaphlan_reference_runner.py,
inside the pinned analysis image. Shell syntax must pass too.

Verify job exit, stderr, root SUCCESS and SHA256SUMS before summarizing.
Endpoint generation is not a performance conclusion. Median/IQR summaries,
paired differences and plots are a subsequent reporting step; no inferential
claims or confidence intervals follow from descriptive sample IQRs.
