#!/usr/bin/env Rscript

script_arg <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(script_arg), "../.."))
tmp <- tempfile("perturbation_response_report.")
dir.create(tmp)
input <- file.path(tmp, "input")
out <- file.path(tmp, "report")
recovery <- file.path(tmp, "recovery")
dir.create(input)
dir.create(recovery)

write_tsv <- function(d, name) {
  write.table(d, file.path(input, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

cohorts <- c("feng", "zeller", "yachida")
populations <- c("community", "independent")
profilers <- c("kraken2_bracken", "metaphlan4")
targets <- c("Bfrag", "Fnuc", "Pmic")
features <- c(targets, paste0("bystander_", 1:8))
operator <- expand.grid(
  holdout_cohort = cohorts, cohort = cohorts, analysis_population = populations,
  profiler = profilers, target_label = targets, feature = features,
  stringsAsFactors = FALSE
)
operator <- operator[operator$cohort != operator$holdout_cohort, ]
operator$feature_role <- ifelse(operator$feature == operator$target_label, "target", "bystander")
operator$eligible_contexts <- 30
operator$distinct_samples <- 10
operator$distinct_doses <- 6
operator$operator_slope <- ifelse(
  operator$feature_role == "target", .90,
  ((match(operator$feature, features) %% 5) - 2) / 40
)
operator$operator_intercept <- 0
operator$operator_r_squared <- ifelse(operator$feature_role == "target", .92, .35)
operator$mean_signed_error <- operator$operator_slope / 10
operator$mean_absolute_error <- abs(operator$mean_signed_error)
operator$detection_transition_rate <- ifelse(operator$feature_role == "target", .02, .08)
write_tsv(operator, "response_operator.tsv")

superposition <- expand.grid(
  cohort = cohorts, profiler = profilers,
  dose_fraction_nominal = c(.0001, .001, .01, .05),
  feature_role = c("target", "bystander"), stringsAsFactors = FALSE
)
superposition$eligible_contexts <- 20
superposition$composition_only_mae <- .02 + 0.2 * superposition$dose_fraction_nominal
superposition$operator_prediction_mae <- .012 + 0.1 * superposition$dose_fraction_nominal
superposition$operator_prediction_rmse <- 1.2 * superposition$operator_prediction_mae
superposition$operator_r_squared <- .78
superposition$calibration_slope <- .94
superposition$improvement_over_composition <-
  (superposition$composition_only_mae - superposition$operator_prediction_mae) /
  superposition$composition_only_mae
write_tsv(superposition, "superposition_summary.tsv")

certificates <- expand.grid(
  holdout_cohort = cohorts, profiler = profilers, feature = features,
  stringsAsFactors = FALSE
)
certificates$training_cohorts <- vapply(certificates$holdout_cohort, function(x) {
  paste(setdiff(cohorts, x), collapse = ",")
}, character(1))
certificates$eligible_contexts <- 240
certificates$eligible_targets <- 9
certificates$eligible_samples <- 40
certificates$baseline_absent_contexts <- 120
certificates$baseline_present_contexts <- 120
offset <- (match(certificates$feature, features) - 1) / 100
certificates$unexpected_appearance_rate <- pmin(.2, .02 + offset)
certificates$unexpected_dropout_rate <- pmin(.2, .03 + offset)
certificates$quantitative_instability_rate <- pmin(.2, .04 + offset)
certificates$mean_absolute_error <- .01 + offset / 10
certificates$hallucination_reliability <- 1 - certificates$unexpected_appearance_rate
certificates$stability_reliability <- 1 - certificates$unexpected_dropout_rate
certificates$quantitative_reliability <- 1 - certificates$quantitative_instability_rate
certificates$overall_reliability <- rowMeans(certificates[c(
  "hallucination_reliability", "stability_reliability", "quantitative_reliability")])
certificates$reliability_ci_lower <- pmax(0, certificates$overall_reliability - .04)
certificates$reliability_ci_upper <- pmin(1, certificates$overall_reliability + .04)
certificates$evaluation_contexts <- 60
certificates$evaluation_overall_reliability <- pmax(
  0, certificates$overall_reliability - rep(c(.01, .03), length.out = nrow(certificates)))
# Component-level NA is scientifically valid when its eligibility denominator is zero.
certificates$unexpected_appearance_rate[1] <- NA_real_
certificates$hallucination_reliability[1] <- NA_real_
certificates$baseline_absent_contexts[1] <- 0
certificates$mean_absolute_error[2] <- NA_real_
write_tsv(certificates, "reliability_certificates.tsv")

detection <- expand.grid(
  cohort = cohorts, analysis_population = populations, profiler = profilers,
  target_label = targets, feature_role = c("target", "bystander"),
  dose_level = 1:3, stringsAsFactors = FALSE
)
detection$holdout_cohort <- "none"
detection$dose_fraction_nominal <- c(.0001, .001, .01)[detection$dose_level]
detection$eligible_contexts <- 20
detection$baseline_absent_contexts <- 10
detection$baseline_present_contexts <- 10
detection$remained_absent <- 9
detection$unexpected_appearances <- 1
detection$remained_detected <- 9
detection$unexpected_dropouts <- 1
detection$appearance_rate <- .10
detection$dropout_rate <- .10
write_tsv(detection, "detection_response_summary.tsv")

quantitative <- detection[c("holdout_cohort", "cohort", "analysis_population", "profiler",
                            "target_label", "feature_role", "dose_level",
                            "dose_fraction_nominal", "eligible_contexts")]
quantitative$mean_expected_abundance <- .02
quantitative$mean_observed_abundance <- .021
quantitative$mean_signed_error <- .001
quantitative$mean_absolute_error <- .002
quantitative$rmse <- .003
quantitative$mean_log2_response_error <- .05
quantitative$quantitative_instability_rate <- .08
write_tsv(quantitative, "quantitative_response_summary.tsv")

panel <- expand.grid(
  holdout_cohort = cohorts, profiler = profilers,
  panel_size = c(1, 2, 3, 5, 10), panel_replicate = 1:4,
  stringsAsFactors = FALSE
)
panel$rank_correlation_full_panel <- pmin(.99, .35 + .065 * panel$panel_size)
panel$heldout_reliability_mae <- pmax(.01, .16 - .012 * panel$panel_size)
panel$high_risk_taxon_recall <- pmin(1, .25 + .075 * panel$panel_size)
write_tsv(panel, "panel_saturation.tsv")

heldout <- expand.grid(
  holdout_cohort = cohorts, profiler = profilers,
  feature_role = c("target", "bystander"), stringsAsFactors = FALSE
)
heldout$training_cohorts <- vapply(heldout$holdout_cohort, function(x) {
  paste(setdiff(cohorts, x), collapse = ",")
}, character(1))
heldout$eligible_contexts <- 120
heldout$distinct_samples <- 20
heldout$composition_only_mae <- .03
heldout$operator_prediction_mae <- .018
heldout$operator_prediction_rmse <- .022
heldout$operator_r_squared <- .72
heldout$improvement_over_composition <- .40
heldout$validation_type <- "leave-one-cohort-out"
write_tsv(heldout, "heldout_operator_validation.tsv")

cohort_holdout <- expand.grid(holdout_cohort = cohorts, profiler = profilers,
                              stringsAsFactors = FALSE)
cohort_holdout$training_cohorts <- vapply(cohort_holdout$holdout_cohort, function(x) {
  paste(setdiff(cohorts, x), collapse = ",")
}, character(1))
cohort_holdout$features_trained <- 100
cohort_holdout$features_evaluated <- 95
cohort_holdout$features_shared <- 90
cohort_holdout$spearman_reliability <- .82
cohort_holdout$mean_absolute_reliability_gap <- .06
cohort_holdout$mean_training_reliability <- .88
cohort_holdout$mean_evaluation_reliability <- .84
write_tsv(cohort_holdout, "cohort_holdout_summary.tsv")

biomarker <- expand.grid(
  holdout_cohort = cohorts, profiler = profilers, feature = paste0("marker_", 1:24),
  stringsAsFactors = FALSE
)
biomarker$overall_reliability <- rep(seq(.05, .98, length.out = 24),
                                    length.out = nrow(biomarker))
biomarker$replication_status <- ifelse(
  biomarker$overall_reliability > .55, "DIRECTIONALLY_REPLICATED", "NOT_SIGNIFICANT_ELSEWHERE")
write_tsv(biomarker, "biomarker_replication.tsv")

target_recovery <- expand.grid(
  cohort = cohorts, condition = c("Control", "CRC"), profiler = profilers,
  target_feature = targets, nominal_target_fraction = c(.0001, .001, .01),
  stringsAsFactors = FALSE
)
target_recovery$samples <- 10
target_recovery$good_fraction <- .8
target_recovery$average_fraction <- .15
target_recovery$poor_fraction <- .05
target_recovery$median_observed_over_expected <- 1
target_recovery$median_absolute_relative_error <- .05
write.table(target_recovery, file.path(recovery, "target_recovery_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
fnuc_prevalence <- expand.grid(cohort = cohorts, condition = c("Control", "CRC"),
                               profiler = profilers, stringsAsFactors = FALSE)
fnuc_prevalence$samples <- 10
fnuc_prevalence$detected_samples <- 7
fnuc_prevalence$baseline_prevalence <- .7
fnuc_prevalence$median_baseline_abundance <- .001
write.table(fnuc_prevalence, file.path(recovery, "fnuc_baseline_prevalence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
fnuc_classes <- expand.grid(cohort = cohorts, condition = c("Control", "CRC"),
  profiler = profilers, nominal_target_fraction = c(.0001, .001, .01),
  recovery_class = c("Good", "Average", "Poor / missed"), stringsAsFactors = FALSE)
fnuc_classes$samples <- 10
fnuc_classes$fraction <- c(.7, .2, .1)
write.table(fnuc_classes, file.path(recovery, "fnuc_recovery_classes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
thresholds <- expand.grid(cohort = cohorts, analysis_population = populations,
  target_label = targets, profiler = profilers, stringsAsFactors = FALSE)
thresholds$minimum_called_fraction <- .001
thresholds$ever_called <- 1
thresholds$tested_doses <- 6
write.table(thresholds, file.path(recovery, "biomarker_detection_thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

script <- file.path(root, "analysis_v2", "scripts", "make_perturbation_response_report.R")
status <- system2("Rscript", c(script, "--input-root", input,
                                "--target-recovery-root", recovery, "--outdir", out,
                                "--report-status", "DEVELOPMENT_ONLY"))
expected_figures <- c(
  "response_operator_crosstalk", "response_operator_performance_map",
  "superposition_error_comparison", "superposition_improvement_by_dose",
  "measurement_reliability_landscape", "measurement_reliability_transport",
  "detection_instability_by_dose", "quantitative_error_by_dose",
  "sentinel_panel_saturation", "heldout_operator_generalization",
  "cohort_holdout_reliability_summary",
  "biomarker_replication_by_reliability", "implanted_target_recovery_map",
  "fnuc_baseline_detection", "fnuc_recovery_classes",
  "paired_biomarker_detection_threshold"
)
stopifnot(
  status == 0,
  file.exists(file.path(out, "SUCCESS")),
  file.exists(file.path(out, "DEVELOPMENT_ONLY.txt")),
  file.exists(file.path(out, "FIGURE_GUIDE.md")),
  file.exists(file.path(out, "provenance", "perturbation_response_report.sha256")),
  all(file.exists(file.path(out, "figures", paste0(expected_figures, ".pdf")))),
  all(file.exists(file.path(out, "figures", paste0(expected_figures, ".png")))),
  all(file.exists(file.path(out, "figure_source", paste0(expected_figures, ".tsv"))))
)

# Core tables are mandatory, optional analyses are not.
minimal_input <- file.path(tmp, "minimal_input")
minimal_out <- file.path(tmp, "minimal_report")
dir.create(minimal_input)
for (name in c("response_operator.tsv", "superposition_summary.tsv",
               "reliability_certificates.tsv")) {
  file.copy(file.path(input, name), file.path(minimal_input, name))
}
minimal_status <- system2("Rscript", c(script, "--input-root", minimal_input,
                                        "--outdir", minimal_out,
                                        "--report-status", "DEVELOPMENT_ONLY"))
stopifnot(
  minimal_status == 0,
  file.exists(file.path(minimal_out, "SUCCESS")),
  nrow(read.delim(file.path(minimal_out, "diagnostics", "analysis_availability.tsv"))) == 9
)
message("[PASS] perturbation-response report fixture")
