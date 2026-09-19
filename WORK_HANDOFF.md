# Project handoff: controlled taxonomic perturbation and CRC biomarkers

Last updated: 2026-09-19. Read this first when resuming the project.

## Status note, 18 September 2026

`analyze_perturbation_response.py` now obeys the selected `--reference-scale` on
every quantitative path: the retained baseline is the supplied selected value
(never reconstructed as `(1-F)o`), the response operator and the cross-cohort
holdout prediction regress on a per-perturbation selected-scale implanted driver
(`q_it/D_i` for MetaPhlAn under `profiler_scale`), and
`exact_superposition_rows` selects one internally consistent column triple and
fails closed on the wrong one. Detection, Bracken and the read-proportional
columns are unchanged. Guarded by 17 hand-computed tests in
`analysis_v2/tests/test_builder_analyzer_profiler_scale.py`.

Codex re-verification passed on 19 September 2026, including the fail-closed
exactly-one-driver guard. The cluster handoff
(`analysis_v2/CLUSTER_HANDOFF_GEFF_PROPAGATION.md`) is ready for **André's
DEVELOPMENT_ONLY execution**. `mgcv` is absent locally, so
`fit_continuous_dose_response.R` has been parse-checked only and has never been
executed; the development run must establish that path before any definitive
claim.

**PENDING, not started:** the Figure 6 implanted-taxon exclusion defect —
community profiles exclude only the focal target rather than all ten implanted
taxa, community profiles may be counted repeatedly across `target_label` rows,
and feature aliases must be canonicalized before target exclusion. Logged as
`OPEN` in `analysis_v2/METHODS_DECISION_LOG.md`; it needs its own task.

## Scientific story

The preprint mainly asked whether artificial read spike-ins could calibrate
taxonomic profilers and remove apparent false-positive biomarkers. The revised
project asks a broader question:

> Controlled read perturbations reveal profiler-specific detection limits,
> quantitative bias, taxonomic cross-talk, mixture non-additivity, and
> biomarker fragility that ordinary case-control benchmarking cannot identify.

Known paired reads are added to already sequenced stool metagenomes. Ground
truth is implanted read identity, pair count, and final-library read fraction;
it is not cellular abundance, biomass, or extraction efficiency. The cohorts
are FengQ 2015, YachidaS 2019, and ZellerG 2014. Display conditions strictly as
`Control`, `Adenoma`, `CRC`. Profilers are Kraken2 + Bracken and MetaPhlAn 4.
Profiler-native abundances are primary; denominator harmonization is a
sensitivity analysis.

## What changed from the preprint

1. Added complete Yachida and moved from two to three cohorts.
2. Replaced primarily unpaired comparisons with paired, baseline-adjusted
   models using biological sample as the replication unit.
3. Separated independent single-taxon and ten-member community perturbations.
4. Added detection and continuous dose-response analyses.
5. Added a target-to-bystander response operator (cross-talk atlas).
6. Added community-superposition tests asking whether independently learned
   responses predict mixtures better than composition-only dilution.
7. Expanded analysis beyond implanted taxa and baseline biomarkers to the full
   reported feature universe.
8. Reframed the old artefact filter: disappearance after a spike does not alone
   prove a false positive because dilution, thresholding, compositional effects,
   and profiler instability are alternatives.
9. Added disease-biomarker retention, effect fidelity, direction stability,
   induced calls, and continuous perturbation-reliability annotations.
10. Kept cross-cohort replication as an independent annotation, not truth.
11. Added calibration-to-biomarker linkage and an ideal counterfactual
    difference-in-differences analysis.
12. Prefer continuous `log2(observed / expected)` recovery. The 10%/50%
    recovery classes are secondary descriptive summaries.
13. Added explicit `DEVELOPMENT_ONLY`/`DEFINITIVE` status, diagnostics,
    manifests, checksums, and fail-closed completion markers.

## Evidence status

- Yachida strict production is complete and sealed: 201 samples, 201 baseline,
  1,407 community, and 1,800 independent profiles.
- The transferred three-cohort downstream run combines complete Yachida with
  historical Feng/Zeller profiles. It is `DEVELOPMENT_ONLY` engineering and
  figure-development evidence.
- Strict-production Feng and Zeller must replace those historical inputs before
  definitive manuscript claims.
- Slurm `COMPLETED` is insufficient: require stage `SUCCESS`, diagnostics,
  receipts, and checksum seals.

### Live upstream production snapshot (2026-09-15)

- **Yachida:** complete and dataset-sealed. The production seal checksum was
  verified successfully. This cohort is ready as upstream evidence.
- **Zeller:** strict production is still running from
  `config/zeller.strict-production.env`. At the latest observed checkpoint,
  138 of 156 samples had both `.verified` markers and retained-output receipts.
  The rolling submission covers all 156 manifest positions. No failure was
  reported in the most recent check of the active final-generation jobs.
- **Feng:** strict-production environment and manifests passed preflight, but
  no current completion count was verified in this local handoff. Treat Feng
  strict production as incomplete/unknown until it is audited on the cluster.

The Zeller rolling ledger is:

`work/crc_production/logs/production_manifest_rolling_20260911T233426Z/jobs.tsv`

Its four generations are download/compute job pairs
`3096721/3096722`, `3096723/3096724`, `3096725/3096774`, and
`3096775/3096776`. Targeted timeout retries used later job IDs and must be
included in the operational audit rather than inferred from this static list.

Never assume this dated snapshot is still current. Recount persistent markers
and query Slurm before taking action.

## Upstream production architecture and recovery

CRC production uses `datasets/submit_crc_rolling.sh`. The intended design is:

- 42 compute lanes;
- no more than 8 simultaneous ENA downloads;
- verified staging before compute;
- per-lane dependencies between generations;
- retained-output verification before disposable input deletion;
- recovery from persistent state rather than resubmitting completed samples.

The generation dependencies proved too conservative when Slurm retained a
dependency after the preceding task in that lane had completed. A safe
operational recovery is to release **download tasks only** when download
capacity is available **and the preceding compute task for the same lane is
confirmed `COMPLETED`**. Never clear the matching compute dependency before
that download succeeds.

Example for pending download tasks:

```bash
download_job=3096775
previous_compute_job=3096774
running_downloads="$(squeue -h -r -j 3096723,3096725,3096775 -t R -o '%i' | wc -l)"
slots=$((8 - running_downloads))

mapfile -t release_tasks < <(
  for lane in $(squeue -h -r -j "$download_job" -t PD -o '%K'); do
    state="$(sacct -X -n -P -j "${previous_compute_job}_${lane}" \
      --format=State | cut -d'|' -f1 | head -1 | sed 's/+.*//')"
    test "$state" = COMPLETED && echo "${download_job}_${lane}"
  done | head -n "$slots"
)

for task in "${release_tasks[@]}"; do
  scontrol update JobId="$task" Dependency=
done
```

This is an operator recovery procedure for one generation pair, not a
replacement for the submitter. Set the two job IDs for the generation being
examined, verify proposed task identities and capacity, and release nothing if
the predecessor state is ambiguous. Corresponding compute jobs remain
dependency-gated.

Current completion must be measured from persistent state:

```bash
source config/zeller.strict-production.env
find "$CRC_STATE_DIR/samples" -maxdepth 1 -type f -name '*.verified' | wc -l
find "$CRC_STATE_DIR/samples" -maxdepth 1 -type f -name '*.retained_outputs.tsv' | wc -l
```

Both counts must reach 156 and agree before the final Zeller cohort audit.

## Taxon identity caveat

Do not collapse profiler taxa by broad prefixes. Mapping
`Fusobacterium nucleatum_E` to `Fusobacterium nucleatum` inflated Kraken2
baseline prevalence and was removed. Corrected provisional Fnuc figures use the
uncollapsed response atlas. A definitive atlas needs a rebuild after the alias
policy is frozen.

Canonical spike order is:

`Bfrag, Csym, Dpne, Fnuc, Hhat, Pmic, Pana, Psto, Porp, Pint`

Use `spikes/spike_panel.tsv` plus the legacy `../spike_taxon_aliases.csv`.
Important aliases include `Allisonella pneumosintes -> Dpne`,
`Hungatella_A hathewayi_A -> Hhat`, and
`Clostridium_Q symbiosum -> Csym`.

## Directory map

The parent `codex_cleaning/` contains transferred evidence and local figure
products. This repository snapshot contains source code and more `work/`
outputs.

### Transferred evidence in the parent directory

- `../three_cohort_analysis_development_bundle_20260914/`: canonical input,
  endpoints, models, reports, and provenance.
- `../three_cohort_complete_figure_suite_dev_20260914/`: consolidated earlier
  three-cohort figure suite.
- `../manuscript_figure_candidates_three_cohort_dev_20260914/`: manuscript
  figure candidates.
- `../three_cohort_disease_robustness_dev_20260914/`: disease transition and
  robustness evidence.
- `../three_cohort_disease_robustness_figures_dev_20260914_fixed2/`: latest
  corrected disease-robustness figures from that cycle.
- `../three_cohort_perturbation_reliability_dev_20260914/`: reliability results.
- `../three_cohort_cross_cohort_artifact_filter_dev_20260914/`: experimental
  filter screen; not validated truth.

### Current atlas and poster work under `work/`

- `work/analysis_v2_perturbation_response_atlas_three_cohort_dev_20260915/`:
  corrected uncollapsed response atlas.
- Its `report.pre_final_superposition/`: analysis tables used for later
  report-only regeneration.
- `work/target_recovery_poster_three_cohort_community_corrected_v3/`: corrected
  community target-recovery summaries.
- Atlas `report_poster_fnuc_community_v7/`: latest verified Fnuc baseline and
  pooled recovery-class report.
- `../tax_biomarker_benchmarking_for_biomarker/work/target_recovery_continuous_figure1_v4/`:
  current continuous Figure 1, with corrected aliases and fixed poster taxon
  order. Its current plotting script is in the same adjacent development tree
  at `../tax_biomarker_benchmarking_for_biomarker/analysis_v2/scripts/`.

These are development products unless explicitly sealed as definitive.

## Reading order for another agent

1. `WORK_HANDOFF.md`
2. `REPRODUCING.md`
3. `analysis_v2/README.md`
4. `analysis_v2/STATISTICAL_ANALYSIS_PLAN.md` and `ANALYSIS_POLICY.tsv`
5. `analysis_v2/THREE_COHORT_DEVELOPMENT.md`
6. `analysis_v2/PERTURBATION_RESPONSE_ATLAS.md`
7. `analysis_v2/BIOMARKER_PROPAGATION.md`, `DISEASE_BIOMARKER_MODEL.md`, and
   `CALIBRATION_AWARE_ROBUSTNESS.md`
8. `analysis_v2/REPORTING.md` and `THREE_COHORT_PUBLICATION_REPORT.md`
9. `analysis_v2/METHODS_DECISION_LOG.md`
10. The appropriate Yachida or CRC definitive runbook before cluster work.

## Rebuild the mixed three-cohort downstream run on the cluster

Use the authoritative cluster checkout:

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

Keep `YACHIDA_MANIFEST` exported. Use a new output directory. Confirm the final
package and every required report have non-empty `SUCCESS` files.

## Rebuild the perturbation-response atlas

Read `analysis_v2/run_perturbation_response_atlas.sh` and
`PERTURBATION_RESPONSE_ATLAS.md` first, then use a validated canonical input:

```bash
export INPUT_ROOT=/path/to/validated/three_cohort_input
export OUTDIR=/new/empty/atlas_output
export REPORT_STATUS=DEVELOPMENT_ONLY
bash analysis_v2/run_perturbation_response_atlas.sh
```

Required variables may evolve; run fixtures first. Definitive execution must
use final Feng/Zeller inputs and the frozen identity policy.

## Regenerate the current Fnuc report without rebuilding the atlas

```bash
Rscript analysis_v2/scripts/make_perturbation_response_report.R \
  --input-root work/analysis_v2_perturbation_response_atlas_three_cohort_dev_20260915/report.pre_final_superposition \
  --target-recovery-root work/target_recovery_poster_three_cohort_community_corrected_v3 \
  --outdir work/new_empty_report_directory \
  --report-status DEVELOPMENT_ONLY
```

Expected poster outputs include `fnuc_baseline_detection`,
`fnuc_recovery_classes`, and `implanted_target_recovery_map` as PNG and PDF.
Baseline detection uses community samples, profiler rows, cohort columns, and
condition colors. Recovery classes pool cohorts into profiler rows and three
low-dose columns.

## Regenerate continuous Figure 1

```bash
Rscript ../tax_biomarker_benchmarking_for_biomarker/analysis_v2/scripts/make_target_recovery_continuous_plot.R \
  --input work/target_recovery_poster_three_cohort_community_corrected_v3/target_recovery_summary.tsv \
  --outdir work/new_empty_continuous_figure_directory
```

If the summary is in an adjacent transferred tree, locate it with
`find .. -name target_recovery_summary.tsv`. The output must contain PNG, PDF,
the source summary, `DEVELOPMENT_ONLY.txt`, and `SUCCESS`.

## Validation checklist

1. Record source commit and analysis-container checksum.
2. Record whether every cohort input is historical or strict production.
3. Validate canonical input and profiler semantics.
4. Require non-empty stage `SUCCESS` files.
5. Inspect diagnostics and figure-source tables, not only PNGs.
6. Check condition and canonical spike ordering.
7. Audit aliases against controlled identities.
8. Preserve `NA` when an estimand is undefined; never replace it with zero.
   Disease retention is undefined when a stratum has no baseline disease calls.
9. Treat recovery thresholds as secondary; retain continuous recovery and
   detection results.
10. Keep development and definitive directories separate and immutable.

## Immediate next work

1. Finalize the poster using corrected three-cohort development figures.
2. Audit continuous Figure 1 source data and caption.
3. Select a minimal poster panel set without conflating recovery, cross-talk,
   and biomarker robustness.
4. Finish and seal strict-production Feng and Zeller.
5. Rebuild all downstream stages from definitive inputs.
6. Re-run taxon-collapse sensitivity with the corrected alias policy.
7. Freeze the manuscript figure set and update the decision log.

When delegating, specify whether a task is **code only**, **development-output
interpretation**, or **definitive execution**. Never let development results be
silently promoted into manuscript evidence.
