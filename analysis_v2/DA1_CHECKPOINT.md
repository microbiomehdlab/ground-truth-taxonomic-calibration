# DA1 eligibility and baseline-family checkpoint

No statistical fits. Shared launcher: run_three_cohort_da1_checkpoint.sbatch.
Requires PROJECT, ANALYSIS_SIF, verified DA_INVENTORY_ROOT and fresh
DA1_CHECKPOINT_ROOT. Reads the existing inventory and rechecks its checksums
and all recorded source hashes; does not change the DA2/DA3 inventory.

Separate Adenoma/Control and CRC/Control comparisons, cohorts and profilers.
Primary complete cases require finite age in (0,120] years and exact Female/Male
sex coding. Invalid values are recorded as exclusions, never imputed. Missing
primary columns stop execution. BMI is retained but does not determine primary
eligibility; BMI sensitivity eligibility remains a separate task. No automatic
batch adjustment. The recorded processing batch is not proof of sequencing batch.

Check rank four and positive residual degrees of freedom for intercept, disease
indicator (Control reference), age and sex (Female reference). Gaussian
elimination checks structural rank; this is not a condition-number, biological
confounding, residual-diagnostic or statistical-calibration assessment.

Freeze baseline prevalence >=10% in each eligible comparison pool, union the ten
exact frozen target aliases, using unchanged native species abundance. Each
profiler keeps its own family and denominator. No post-spike selection/reclosure.

Outputs per cohort: eligibility.tsv (including exclusions and raw age/sex/BMI),
model_design.tsv, feature_families.tsv, input_hashes.tsv, SUCCESS and SHA256SUMS.
Root records commit/image hash and recursive checksums. Primary metadata is
reviewed before modeling; checkpoint PASS is not definitive model readiness.

Cluster metadata inspection on 2026-10-01 found complete age/sex for all 511
samples, Female/Male coding, and five missing Zeller BMI values. DA2/DA3 families
contain 10 targets each, with 3214–3574 Bracken features and 474–795 MetaPhlAn
features. These observations are user-provided cluster output, not local data.
DA1 family counts require the new checkpoint and cannot be inferred from these.

Tests: test_da1_checkpoint.py covers missing/nonfinite/out-of-range age, invalid
sex, full versus confounded rank, integrated verified input/family construction,
existing-output refusal and changed-source refusal before writing.
