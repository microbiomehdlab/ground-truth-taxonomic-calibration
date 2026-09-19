#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}

input_root <- value("--input-root")
outdir <- value("--outdir")
recovery_root <- value("--target-recovery-root")
report_status <- value("--report-status", "DEVELOPMENT_ONLY")
if (is.null(input_root) || is.null(outdir)) {
  stop("Required arguments: --input-root ROOT --outdir OUT [--report-status DEVELOPMENT_ONLY|DEFINITIVE]")
}
if (!report_status %in% c("DEVELOPMENT_ONLY", "DEFINITIVE")) {
  stop("--report-status must be DEVELOPMENT_ONLY or DEFINITIVE")
}
if (!dir.exists(input_root)) stop("Input root does not exist: ", input_root)
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE))) {
  stop("OUTDIR must be new or empty")
}

locate_table <- function(name, required = FALSE) {
  candidates <- c(file.path(input_root, name), file.path(input_root, "tables", name))
  found <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
  if (!length(found)) {
    if (required) stop("Missing required response table: ", name)
    return(NULL)
  }
  normalizePath(found[[1L]])
}

read_table <- function(name, required = FALSE) {
  path <- locate_table(name, required)
  if (is.null(path)) return(NULL)
  ans <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = c("NA", "NaN", ""))
  if (!nrow(ans)) stop("Empty response table: ", path)
  attr(ans, "source_path") <- path
  ans
}

read_recovery_table <- function(name) {
  if (is.null(recovery_root)) return(NULL)
  path <- file.path(recovery_root, name)
  if (!file.exists(path) || file.info(path)$size <= 0) return(NULL)
  ans <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = c("NA", "NaN", ""))
  if (!nrow(ans)) stop("Empty target-recovery table: ", path)
  attr(ans, "source_path") <- normalizePath(path)
  ans
}

require_columns <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

numeric_column <- function(d, column, label, lower = -Inf, upper = Inf,
                           allow_na = FALSE) {
  x <- suppressWarnings(as.numeric(d[[column]]))
  bad <- !is.finite(x) | x < lower | x > upper
  if (allow_na) bad <- !is.na(x) & bad
  if (any(bad)) stop(label, " contains invalid ", column)
  x
}

write_tsv <- function(d, path) {
  write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_rms <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sqrt(sum(x[keep]^2 * w[keep]) / sum(w[keep]))
}

median_iqr <- function(x) {
  data.frame(y = median(x, na.rm = TRUE),
             ymin = unname(quantile(x, .25, na.rm = TRUE)),
             ymax = unname(quantile(x, .75, na.rm = TRUE)))
}

group_apply <- function(d, keys, fun) {
  encoded <- do.call(paste, c(lapply(d[keys], function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "NA"
    x
  }), sep = "\r"))
  pieces <- split(seq_len(nrow(d)), encoded, drop = TRUE)
  rows <- lapply(pieces, function(index) {
    key_values <- d[index[[1L]], keys, drop = FALSE]
    cbind(key_values, fun(d[index, , drop = FALSE]))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

profiler_names <- c(
  kraken2_bracken = "Kraken2 + Bracken",
  metaphlan4 = "MetaPhlAn 4"
)
decorate <- function(d) {
  d$Profiler <- unname(profiler_names[d$profiler])
  d$Profiler[is.na(d$Profiler)] <- d$profiler[is.na(d$Profiler)]
  if ("cohort" %in% names(d)) d$Cohort <- tools::toTitleCase(d$cohort)
  if ("holdout_cohort" %in% names(d)) {
    d$Holdout <- tools::toTitleCase(d$holdout_cohort)
    d$Holdout[is.na(d$Holdout) | !nzchar(d$Holdout)] <- "None"
  }
  if ("analysis_population" %in% names(d)) {
    d$Population <- tools::toTitleCase(d$analysis_population)
  }
  if ("condition" %in% names(d)) {
    d$Condition <- factor(d$condition,
      levels = c("Control", "Adenoma", "CRC", "colorectal carcinoma"))
  }
  d
}

short_feature <- function(x, width = 38L) {
  ifelse(nchar(x) > width, paste0(substr(x, 1L, width - 1L), "\u2026"), x)
}

for (directory in c("figures", "figure_source", "tables", "diagnostics", "provenance")) {
  dir.create(file.path(outdir, directory), recursive = TRUE, showWarnings = FALSE)
}

operator <- read_table("response_operator.tsv", required = TRUE)
superposition <- read_table("superposition_summary.tsv", required = TRUE)
certificates <- read_table("reliability_certificates.tsv", required = TRUE)
detection <- read_table("detection_response_summary.tsv")
quantitative <- read_table("quantitative_response_summary.tsv")
heldout <- read_table("heldout_validation.tsv")
if (is.null(heldout)) heldout <- read_table("heldout_operator_validation.tsv")
cohort_holdout <- read_table("cohort_holdout_summary.tsv")
panel <- read_table("panel_saturation.tsv")
biomarker <- read_table("biomarker_replication.tsv")
target_recovery <- read_recovery_table("target_recovery_summary.tsv")
fnuc_prevalence <- read_recovery_table("fnuc_baseline_prevalence.tsv")
fnuc_classes <- read_recovery_table("fnuc_recovery_classes.tsv")
biomarker_thresholds <- read_recovery_table("biomarker_detection_thresholds.tsv")
input_tables <- Filter(Negate(is.null), list(operator, superposition, certificates,
                                              detection, quantitative, heldout, cohort_holdout,
                                              panel, biomarker, target_recovery,
                                              fnuc_prevalence, fnuc_classes,
                                              biomarker_thresholds))
input_paths <- unname(vapply(input_tables, attr, character(1), which = "source_path"))

require_columns(operator, c(
  "holdout_cohort", "cohort", "analysis_population", "profiler", "target_label",
  "feature", "feature_role", "eligible_contexts", "operator_slope",
  "operator_intercept", "operator_r_squared", "mean_signed_error",
  "mean_absolute_error", "detection_transition_rate"
), "response_operator.tsv")
operator$eligible_contexts <- numeric_column(operator, "eligible_contexts", "response_operator.tsv", 1)
for (field in c("operator_slope", "operator_intercept", "mean_signed_error",
                "mean_absolute_error", "detection_transition_rate")) {
  operator[[field]] <- numeric_column(operator, field, "response_operator.tsv")
}
operator$operator_r_squared <- numeric_column(operator, "operator_r_squared",
                                               "response_operator.tsv", allow_na = TRUE)
if (any(operator$mean_absolute_error < 0) ||
    any(operator$detection_transition_rate < 0 | operator$detection_transition_rate > 1)) {
  stop("response_operator.tsv contains impossible error or transition values")
}
operator <- decorate(operator)
operator$is_target <- grepl("target|implant", tolower(operator$feature_role)) &
  !grepl("bystander|off", tolower(operator$feature_role))

require_columns(superposition, c(
  "cohort", "profiler", "dose_fraction_nominal", "feature_role",
  "eligible_contexts", "composition_only_mae", "operator_prediction_mae",
  "operator_prediction_rmse", "operator_r_squared", "calibration_slope",
  "improvement_over_composition"
), "superposition_summary.tsv")
if (!"analysis_population" %in% names(superposition)) {
  # Superposition is intrinsically evaluated on community-mixture profiles.
  superposition$analysis_population <- "community"
}
for (field in c("dose_fraction_nominal", "eligible_contexts", "composition_only_mae",
                "operator_prediction_mae", "operator_prediction_rmse", "operator_r_squared",
                "calibration_slope", "improvement_over_composition")) {
  superposition[[field]] <- numeric_column(superposition, field, "superposition_summary.tsv",
                                           allow_na = field %in% c("operator_r_squared", "calibration_slope"))
}
if (any(superposition$dose_fraction_nominal < 0) ||
    any(superposition$eligible_contexts < 1) ||
    any(superposition$composition_only_mae < 0) ||
    any(superposition$operator_prediction_mae < 0) ||
    any(superposition$operator_prediction_rmse < 0)) {
  stop("superposition_summary.tsv contains impossible dose, count, or error values")
}
superposition <- decorate(superposition)

require_columns(certificates, c(
  "holdout_cohort", "training_cohorts", "profiler", "feature", "eligible_contexts",
  "eligible_targets", "eligible_samples", "unexpected_appearance_rate",
  "unexpected_dropout_rate", "quantitative_instability_rate", "mean_absolute_error",
  "hallucination_reliability", "stability_reliability", "quantitative_reliability",
  "overall_reliability", "reliability_ci_lower", "reliability_ci_upper",
  "evaluation_contexts", "evaluation_overall_reliability"
), "reliability_certificates.tsv")
for (field in c("eligible_contexts", "eligible_targets", "eligible_samples", "evaluation_contexts")) {
  certificates[[field]] <- numeric_column(
    certificates, field, "reliability_certificates.tsv", 0,
    allow_na = field == "evaluation_contexts"
  )
}
component_fields <- c("unexpected_appearance_rate", "unexpected_dropout_rate",
                      "quantitative_instability_rate", "hallucination_reliability",
                      "stability_reliability", "quantitative_reliability")
for (field in component_fields) {
  certificates[[field]] <- numeric_column(certificates, field,
                                          "reliability_certificates.tsv", 0, 1,
                                          allow_na = TRUE)
}
for (field in c("overall_reliability", "reliability_ci_lower", "reliability_ci_upper")) {
  certificates[[field]] <- numeric_column(certificates, field,
                                          "reliability_certificates.tsv", 0, 1)
}
certificates$evaluation_overall_reliability <- numeric_column(
  certificates, "evaluation_overall_reliability", "reliability_certificates.tsv",
  0, 1, allow_na = TRUE)
certificates$mean_absolute_error <- numeric_column(certificates, "mean_absolute_error",
                                                   "reliability_certificates.tsv", 0,
                                                   allow_na = TRUE)
if (any(certificates$reliability_ci_lower > certificates$overall_reliability |
        certificates$reliability_ci_upper < certificates$overall_reliability)) {
  stop("Reliability confidence interval does not contain its estimate")
}
certificates <- decorate(certificates)

theme_set(theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
  legend.position = "bottom", plot.title.position = "plot"
))

plots <- list()
figure_sources <- list()

# Response operator: retain the strongest recurrent bystander edges per profiler.
bystander <- operator[!operator$is_target & is.finite(operator$operator_slope), , drop = FALSE]
if (!nrow(bystander)) stop("The response operator contains no bystander rows")
operator_atlas <- group_apply(
  bystander, c("Profiler", "target_label", "feature"),
  function(x) data.frame(
    operator_slope = weighted_mean(x$operator_slope, x$eligible_contexts),
    eligible_contexts = sum(x$eligible_contexts),
    cohorts = length(unique(x$cohort)),
    stringsAsFactors = FALSE
  )
)
keep_features <- unlist(lapply(split(operator_atlas, operator_atlas$Profiler), function(x) {
  strength <- tapply(abs(x$operator_slope), x$feature, max, na.rm = TRUE)
  names(sort(strength, decreasing = TRUE))[seq_len(min(15L, length(strength)))]
}), use.names = FALSE)
operator_atlas <- operator_atlas[operator_atlas$feature %in% keep_features, , drop = FALSE]
operator_atlas$Feature <- short_feature(operator_atlas$feature)
limit <- max(abs(operator_atlas$operator_slope), na.rm = TRUE)
if (!is.finite(limit) || limit == 0) limit <- 1
plots$response_operator_crosstalk <- ggplot(
  operator_atlas, aes(target_label, Feature, fill = operator_slope)
) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~Profiler, scales = "free_y") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-limit, limit), oob = squish) +
  labs(
    title = "Controlled perturbations identify profiler-specific taxonomic cross-talk",
    subtitle = "Strongest target-to-bystander response slopes; expected compositional change removed",
    x = "Implanted taxon", y = "Reported bystander taxon", fill = "Response\nslope"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
figure_sources$response_operator_crosstalk <- operator_atlas

target_rows <- operator[operator$is_target, , drop = FALSE]
if (!nrow(target_rows)) stop("The response operator contains no implanted-target rows")
performance_keys <- c("holdout_cohort", "cohort", "analysis_population", "profiler",
                      "target_label", "Holdout", "Cohort", "Population", "Profiler")
target_performance <- group_apply(target_rows, performance_keys, function(x) data.frame(
  target_gain_error = abs(weighted_mean(x$operator_slope, x$eligible_contexts) - 1),
  target_contexts = sum(x$eligible_contexts), stringsAsFactors = FALSE
))
bystander_performance <- group_apply(bystander, performance_keys, function(x) data.frame(
  bystander_crosstalk = weighted_rms(x$operator_slope, x$eligible_contexts),
  transition_rate = weighted_mean(x$detection_transition_rate, x$eligible_contexts),
  bystander_contexts = sum(x$eligible_contexts), stringsAsFactors = FALSE
))
performance <- merge(target_performance, bystander_performance, by = performance_keys)
if (!nrow(performance)) stop("Target and bystander operator strata do not overlap")
plots$response_operator_performance_map <- ggplot(
  performance, aes(target_gain_error, bystander_crosstalk, colour = Profiler,
                   shape = Population, size = pmax(1, bystander_contexts))
) +
  geom_point(alpha = 0.76) +
  facet_wrap(~Holdout) +
  scale_colour_manual(values = c("Kraken2 + Bracken" = "#D55E00", "MetaPhlAn 4" = "#0072B2")) +
  scale_size_area(max_size = 7, guide = "none") +
  labs(
    title = "Profiler performance has separable target and bystander dimensions",
    subtitle = "Lower-left is ideal; each point is one held-out cohort, population and target",
    x = "Absolute target-gain error |slope - 1|",
    y = "Bystander cross-talk energy (RMS slope)", colour = "Profiler", shape = "Population"
  )
figure_sources$response_operator_performance_map <- performance

# Independent-to-community superposition.
plots$superposition_error_comparison <- ggplot(
  superposition, aes(composition_only_mae, operator_prediction_mae,
                     colour = Cohort, shape = feature_role,
                     size = pmax(1, eligible_contexts))
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
  geom_point(alpha = 0.78) +
  # Stack profilers so each comparison receives the full figure width.  A
  # fixed coordinate ratio made these panels extremely shallow because the
  # composition-only and operator-error ranges differ substantially.
  facet_grid(Profiler ~ Population) +
  scale_size_area(max_size = 7, guide = "none") +
  labs(
    title = "Single-taxon responses are tested as predictors of community perturbations",
    subtitle = "Points below the diagonal improve on composition-only prediction",
    x = "Composition-only mean absolute error",
    y = "Response-operator mean absolute error", colour = "Cohort", shape = "Feature role"
  )
figure_sources$superposition_error_comparison <- superposition

plots$superposition_improvement_by_dose <- ggplot(
  superposition, aes(100 * dose_fraction_nominal, improvement_over_composition,
                     colour = Cohort,
                     group = interaction(Cohort, feature_role))
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
  stat_summary(fun = median, geom = "line", linewidth = 0.75) +
  stat_summary(fun = median, geom = "point", size = 2) +
  facet_grid(Population ~ Profiler) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Superposition gains reveal when profiler behaviour is approximately additive",
    subtitle = "Positive values indicate lower error than the composition-only expectation",
    x = "Community implanted fraction (%)", y = "Relative error improvement",
    colour = "Cohort"
  )
figure_sources$superposition_improvement_by_dose <- superposition

# Feature-level, target-agnostic measurement certificate.
certificate_landscape <- certificates[
  is.finite(certificates$hallucination_reliability) &
  is.finite(certificates$stability_reliability) &
  is.finite(certificates$quantitative_reliability), , drop = FALSE]
if (!nrow(certificate_landscape)) {
  stop("No reliability certificate has all three components estimable")
}
plots$measurement_reliability_landscape <- ggplot(
  certificate_landscape, aes(hallucination_reliability, stability_reliability,
                    colour = quantitative_reliability,
                    size = pmax(1, evaluation_contexts))
) +
  geom_point(alpha = 0.62) +
  facet_grid(Holdout ~ Profiler) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_colour_viridis_c(limits = c(0, 1), labels = percent_format(), option = "C") +
  scale_size_area(max_size = 6, guide = "none") +
  labs(
    title = "Perturbations provide target-agnostic measurement certificates",
    subtitle = "Features with all three components estimable; evaluated only while a bystander",
    x = "Resistance to unexpected appearance", y = "Resistance to dropout",
    colour = "Quantitative\nreliability"
  ) +
  theme(legend.position = "right")
figure_sources$measurement_reliability_landscape <- certificate_landscape

certificate_summary <- group_apply(
  certificates, c("Holdout", "Profiler"), function(x) data.frame(
    features = nrow(x),
    median_reliability = median(x$overall_reliability, na.rm = TRUE),
    median_evaluation_reliability = median(x$evaluation_overall_reliability, na.rm = TRUE),
    evaluable_features = sum(is.finite(x$evaluation_overall_reliability)),
    stringsAsFactors = FALSE
  )
)
plots$measurement_reliability_transport <- ggplot(
  certificates[is.finite(certificates$evaluation_overall_reliability), , drop = FALSE],
  aes(overall_reliability, evaluation_overall_reliability, colour = Profiler,
      size = pmax(1, evaluation_contexts))
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
  geom_point(alpha = 0.60) +
  facet_wrap(~Holdout) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_colour_manual(values = c("Kraken2 + Bracken" = "#D55E00", "MetaPhlAn 4" = "#0072B2")) +
  scale_size_area(max_size = 6, guide = "none") +
  labs(
    title = "Reliability learned outside a cohort is evaluated in the held-out cohort",
    subtitle = "Agreement with the diagonal indicates transportable feature reliability",
    x = "Training-cohort reliability", y = "Held-out reliability", colour = "Profiler"
  )
figure_sources$measurement_reliability_transport <- certificates
figure_sources$measurement_reliability_summary <- certificate_summary
certificate_estimability <- data.frame(
  metric = c("certificates", "hallucination_component_estimable",
             "stability_component_estimable", "quantitative_component_estimable",
             "all_three_components_estimable", "heldout_reliability_estimable"),
  features = c(nrow(certificates),
    sum(is.finite(certificates$hallucination_reliability)),
    sum(is.finite(certificates$stability_reliability)),
    sum(is.finite(certificates$quantitative_reliability)),
    nrow(certificate_landscape),
    sum(is.finite(certificates$evaluation_overall_reliability))),
  stringsAsFactors = FALSE
)
figure_sources$measurement_reliability_estimability <- certificate_estimability

# Detection and abundance-response diagnostics are emitted when available.
if (!is.null(detection)) {
  require_columns(detection, c("cohort", "analysis_population", "profiler", "feature_role",
    "dose_fraction_nominal", "eligible_contexts", "appearance_rate", "dropout_rate"),
    "detection_response_summary.tsv")
  for (field in c("dose_fraction_nominal", "eligible_contexts", "appearance_rate", "dropout_rate")) {
    detection[[field]] <- numeric_column(detection, field, "detection_response_summary.tsv",
                                         if (grepl("rate", field)) 0 else -Inf,
                                         if (grepl("rate", field)) 1 else Inf,
                                         allow_na = grepl("rate", field))
  }
  detection <- decorate(detection)
  detection_long <- rbind(
    transform(detection, transition = "Unexpected appearance", rate = appearance_rate),
    transform(detection, transition = "Unexpected dropout", rate = dropout_rate)
  )
  plots$detection_instability_by_dose <- ggplot(
    detection_long, aes(100 * dose_fraction_nominal, rate, colour = transition,
                        group = interaction(Cohort, transition, feature_role))
  ) +
    stat_summary(fun = median, geom = "line") +
    stat_summary(fun = median, geom = "point") +
    facet_grid(Population ~ interaction(Cohort, Profiler, sep = "\n")) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
    labs(title = "Detection instability is separated into appearance and dropout",
         x = "Implanted fraction (%)", y = "Median transition probability",
         colour = "Transition")
  figure_sources$detection_instability_by_dose <- detection_long
}

if (!is.null(quantitative)) {
  require_columns(quantitative, c("cohort", "analysis_population", "profiler", "feature_role",
    "dose_fraction_nominal", "eligible_contexts", "mean_absolute_error",
    "quantitative_instability_rate"), "quantitative_response_summary.tsv")
  for (field in c("dose_fraction_nominal", "eligible_contexts", "mean_absolute_error",
                  "quantitative_instability_rate")) {
    quantitative[[field]] <- numeric_column(quantitative, field,
      "quantitative_response_summary.tsv", 0,
      if (field == "quantitative_instability_rate") 1 else Inf,
      allow_na = field == "quantitative_instability_rate")
  }
  quantitative <- decorate(quantitative)
  plots$quantitative_error_by_dose <- ggplot(
    quantitative, aes(100 * dose_fraction_nominal, mean_absolute_error,
                      colour = feature_role,
                      group = interaction(Cohort, feature_role))
  ) +
    stat_summary(fun = median, geom = "line") +
    stat_summary(fun = median, geom = "point") +
    facet_grid(Population ~ interaction(Cohort, Profiler, sep = "\n")) +
    labs(title = "Quantitative error is measured against the exact mixing counterfactual",
         x = "Implanted fraction (%)", y = "Median absolute bounded error",
         colour = "Feature role")
  figure_sources$quantitative_error_by_dose <- quantitative
}

# Optional reduced-panel saturation analysis.
if (!is.null(panel)) {
  require_columns(panel, c("holdout_cohort", "profiler", "panel_size", "panel_replicate",
    "rank_correlation_full_panel", "heldout_reliability_mae", "high_risk_taxon_recall"),
    "panel_saturation.tsv")
  panel$panel_size <- numeric_column(panel, "panel_size", "panel_saturation.tsv", 1)
  for (field in c("rank_correlation_full_panel", "high_risk_taxon_recall")) {
    panel[[field]] <- numeric_column(panel, field, "panel_saturation.tsv", -1, 1)
  }
  panel$heldout_reliability_mae <- numeric_column(panel, "heldout_reliability_mae",
                                                  "panel_saturation.tsv", 0)
  panel <- decorate(panel)
  panel_long <- rbind(
    data.frame(panel[c("Holdout", "Profiler", "panel_size", "panel_replicate")],
               metric = "Rank correlation with full panel",
               value = panel$rank_correlation_full_panel),
    data.frame(panel[c("Holdout", "Profiler", "panel_size", "panel_replicate")],
               metric = "High-risk-taxon recall", value = panel$high_risk_taxon_recall)
  )
  plots$sentinel_panel_saturation <- ggplot(
    panel_long, aes(panel_size, value, colour = Profiler,
                    group = interaction(Profiler, Holdout))
  ) +
    stat_summary(fun = median, geom = "line") +
    stat_summary(fun = median, geom = "point") +
    stat_summary(fun.data = median_iqr, geom = "errorbar", width = 0.15, alpha = 0.5) +
    facet_grid(metric ~ Holdout) +
    scale_y_continuous(labels = percent_format()) +
    scale_x_continuous(breaks = sort(unique(panel$panel_size))) +
    scale_colour_manual(values = c("Kraken2 + Bracken" = "#D55E00", "MetaPhlAn 4" = "#0072B2")) +
    labs(title = "Panel-size saturation tests whether a small sentinel set is sufficient",
         x = "Number of independently implanted taxa", y = "Held-out performance",
         colour = "Profiler")
  figure_sources$sentinel_panel_saturation <- panel_long
}

# Optional generic held-out predictions. Both observation-level and summary schemas are accepted.
if (!is.null(heldout)) {
  require_columns(heldout, c("holdout_cohort", "profiler", "validation_type"),
                  "heldout_validation.tsv")
  heldout <- decorate(heldout)
  if (all(c("observed_reliability", "predicted_reliability") %in% names(heldout))) {
    heldout$observed_reliability <- numeric_column(heldout, "observed_reliability",
                                                   "heldout_validation.tsv", 0, 1)
    heldout$predicted_reliability <- numeric_column(heldout, "predicted_reliability",
                                                    "heldout_validation.tsv", 0, 1)
    plots$heldout_reliability_prediction <- ggplot(
      heldout, aes(observed_reliability, predicted_reliability, colour = Profiler)
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
      geom_point(alpha = 0.60) +
      facet_grid(validation_type ~ Holdout) +
      scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
      scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
      labs(title = "Held-out species and cohorts test target-agnostic generalization",
           x = "Observed reliability", y = "Predicted reliability", colour = "Profiler")
    figure_sources$heldout_reliability_prediction <- heldout
  } else if (all(c("prediction_mae", "rank_correlation") %in% names(heldout))) {
    heldout$prediction_mae <- numeric_column(heldout, "prediction_mae",
                                             "heldout_validation.tsv", 0)
    heldout$rank_correlation <- numeric_column(heldout, "rank_correlation",
                                               "heldout_validation.tsv", -1, 1)
    heldout_long <- rbind(
      data.frame(heldout[c("Holdout", "Profiler", "validation_type")],
                 metric = "Prediction MAE", value = heldout$prediction_mae),
      data.frame(heldout[c("Holdout", "Profiler", "validation_type")],
                 metric = "Rank correlation", value = heldout$rank_correlation)
    )
    plots$heldout_reliability_prediction <- ggplot(
      heldout_long, aes(Holdout, value, colour = Profiler, group = Profiler)
    ) +
      geom_point(position = position_dodge(width = 0.25), size = 2) +
      facet_grid(metric ~ validation_type, scales = "free_y") +
      labs(title = "Held-out species and cohorts test target-agnostic generalization",
           x = "Held-out cohort", y = "Validation metric", colour = "Profiler")
    figure_sources$heldout_reliability_prediction <- heldout_long
  } else if (all(c("composition_only_mae", "operator_prediction_mae",
                   "feature_role", "eligible_contexts") %in% names(heldout))) {
    for (field in c("composition_only_mae", "operator_prediction_mae", "eligible_contexts")) {
      heldout[[field]] <- numeric_column(heldout, field, "heldout_operator_validation.tsv",
                                         0)
    }
    plots$heldout_operator_generalization <- ggplot(
      heldout, aes(composition_only_mae, operator_prediction_mae,
                   colour = Profiler, shape = feature_role,
                   size = pmax(1, eligible_contexts))
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
      geom_point(alpha = 0.75) +
      facet_wrap(~Holdout) +
      scale_size_area(max_size = 7, guide = "none") +
      labs(title = "Training-cohort response operators are evaluated in untouched cohorts",
           subtitle = "Points below the diagonal outperform composition-only prediction",
           x = "Composition-only mean absolute error",
           y = "Held-out operator mean absolute error", colour = "Profiler",
           shape = "Feature role")
    figure_sources$heldout_operator_generalization <- heldout
  } else {
    stop("Held-out validation table has no supported prediction schema")
  }
}

if (!is.null(cohort_holdout)) {
  require_columns(cohort_holdout, c("holdout_cohort", "training_cohorts", "profiler",
    "features_trained", "features_evaluated", "features_shared", "spearman_reliability",
    "mean_absolute_reliability_gap", "mean_training_reliability",
    "mean_evaluation_reliability"), "cohort_holdout_summary.tsv")
  cohort_holdout$spearman_reliability <- numeric_column(
    cohort_holdout, "spearman_reliability", "cohort_holdout_summary.tsv", -1, 1,
    allow_na = TRUE)
  cohort_holdout$mean_absolute_reliability_gap <- numeric_column(
    cohort_holdout, "mean_absolute_reliability_gap", "cohort_holdout_summary.tsv", 0,
    allow_na = TRUE)
  cohort_holdout <- decorate(cohort_holdout)
  cohort_holdout_long <- rbind(
    data.frame(cohort_holdout[c("Holdout", "Profiler", "features_shared")],
               metric = "Spearman reliability", value = cohort_holdout$spearman_reliability),
    data.frame(cohort_holdout[c("Holdout", "Profiler", "features_shared")],
               metric = "Mean absolute reliability gap",
               value = cohort_holdout$mean_absolute_reliability_gap)
  )
  plots$cohort_holdout_reliability_summary <- ggplot(
    cohort_holdout_long,
    aes(Holdout, value, colour = Profiler, group = Profiler,
        size = pmax(1, features_shared))
  ) +
    geom_point(position = position_dodge(width = .3), alpha = .85) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    scale_size_area(max_size = 6, guide = "none") +
    labs(title = "Feature reliability is audited by leave-one-cohort-out transfer",
         subtitle = "Correlation should be high and absolute reliability gap low",
         x = "Held-out cohort", y = "Validation metric", colour = "Profiler")
  figure_sources$cohort_holdout_reliability_summary <- cohort_holdout_long
}

# Optional external biomarker replication analysis.
if (!is.null(biomarker)) {
  require_columns(biomarker, c("profiler", "feature"), "biomarker_replication.tsv")
  reliability_field <- intersect(c("overall_reliability", "reliability_score"), names(biomarker))
  if (!length(reliability_field)) stop("biomarker_replication.tsv lacks a reliability score")
  biomarker$reliability <- numeric_column(biomarker, reliability_field[[1L]],
                                          "biomarker_replication.tsv", 0, 1)
  if ("replicated" %in% names(biomarker)) {
    raw <- tolower(as.character(biomarker$replicated))
    biomarker$replicated_binary <- as.integer(raw %in% c("1", "true", "yes", "replicated"))
  } else if ("replication_status" %in% names(biomarker)) {
    biomarker$replicated_binary <- as.integer(grepl("directionally_replicated|^replicated$",
                                                     tolower(biomarker$replication_status)))
  } else {
    stop("biomarker_replication.tsv lacks replicated or replication_status")
  }
  biomarker <- decorate(biomarker)
  biomarker$Reliability_bin <- cut(biomarker$reliability, breaks = c(0, .25, .5, .75, 1),
                                   include.lowest = TRUE,
                                   labels = c("0-25%", "25-50%", "50-75%", "75-100%"))
  keys <- c("Profiler", "Reliability_bin")
  if ("holdout_cohort" %in% names(biomarker)) keys <- c("Holdout", keys)
  replication_summary <- group_apply(biomarker[!is.na(biomarker$Reliability_bin), ], keys,
    function(x) {
      successes <- sum(x$replicated_binary)
      n <- nrow(x)
      p <- successes / n
      z <- qnorm(.975)
      centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
      half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
      data.frame(biomarkers = n, replicated = successes, replication_rate = p,
                 ci_lower = max(0, centre - half), ci_upper = min(1, centre + half))
    })
  facet_formula <- if ("Holdout" %in% names(replication_summary)) ~Holdout else ~Profiler
  plots$biomarker_replication_by_reliability <- ggplot(
    replication_summary,
    aes(Reliability_bin, replication_rate, colour = Profiler, group = Profiler)
  ) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.12,
                  position = position_dodge(width = 0.25)) +
    geom_point(position = position_dodge(width = 0.25), size = 2.2) +
    facet_wrap(facet_formula) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
    labs(title = "Does perturbation reliability predict external biomarker replication?",
         subtitle = "Wilson intervals; association is evaluated without using the validation cohort",
         x = "Perturbation-reliability tier", y = "Directionally replicated biomarkers",
         colour = "Profiler")
  figure_sources$biomarker_replication_by_reliability <- replication_summary
}

# Optional direct-target recovery and paired biomarker-discovery panels.  These
# are poster-friendly descendants of the original manuscript Figures 1--3.
if (!is.null(target_recovery)) {
  require_columns(target_recovery, c("cohort", "profiler", "target_feature",
    "nominal_target_fraction", "samples", "good_fraction"),
    "target_recovery_summary.tsv")
  target_recovery$nominal_target_fraction <- numeric_column(
    target_recovery, "nominal_target_fraction", "target_recovery_summary.tsv", 0)
  target_recovery$samples <- numeric_column(
    target_recovery, "samples", "target_recovery_summary.tsv", 1)
  target_recovery$good_fraction <- numeric_column(
    target_recovery, "good_fraction", "target_recovery_summary.tsv", 0, 1)
  recovery_overview <- group_apply(target_recovery,
    c("cohort", "profiler", "target_feature", "nominal_target_fraction"),
    function(x) data.frame(samples = sum(x$samples),
      good_fraction = weighted_mean(x$good_fraction, x$samples)))
  recovery_overview <- decorate(recovery_overview)
  recovery_overview$Dose <- factor(
    percent(recovery_overview$nominal_target_fraction, accuracy = .001),
    levels = percent(sort(unique(recovery_overview$nominal_target_fraction)), accuracy = .001))
  recovery_overview$Target <- short_feature(recovery_overview$target_feature, 28L)
  plots$implanted_target_recovery_map <- ggplot(
    recovery_overview, aes(Dose, Target, fill = good_fraction)
  ) +
    geom_tile(colour = "white", linewidth = .25) +
    facet_grid(Profiler ~ Cohort, scales = "free_y", space = "free_y") +
    scale_fill_viridis_c(limits = c(0, 1), labels = percent_format(), option = "C") +
    labs(title = "Implanted-target recovery is dose-, cohort-, and profiler-dependent",
         subtitle = "Fraction recovered within 10% of the exact implanted fraction; conditions pooled descriptively",
         x = "Implanted target fraction", y = "Implanted taxon",
         fill = "Good recovery") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  figure_sources$implanted_target_recovery_map <- recovery_overview
}

if (!is.null(fnuc_prevalence)) {
  require_columns(fnuc_prevalence, c("cohort", "analysis_population", "condition", "profiler", "samples",
    "detected_samples", "baseline_prevalence"), "fnuc_baseline_prevalence.tsv")
  fnuc_prevalence$baseline_prevalence <- numeric_column(
    fnuc_prevalence, "baseline_prevalence", "fnuc_baseline_prevalence.tsv", 0, 1)
  fnuc_prevalence <- decorate(fnuc_prevalence)
  fnuc_prevalence <- fnuc_prevalence[fnuc_prevalence$Population == "Community", , drop = FALSE]
  plots$fnuc_baseline_detection <- ggplot(
    fnuc_prevalence, aes(Condition, baseline_prevalence, fill = Condition)
  ) +
    geom_col(position = position_dodge(width = .75), width = .68) +
    geom_text(aes(label = paste0(detected_samples, "/", samples)),
              position = position_dodge(width = .75), vjust = -.25, size = 3) +
    facet_grid(Profiler ~ Cohort) +
    scale_y_continuous(limits = c(0, 1.08), labels = percent_format()) +
    scale_fill_manual(values = c(Control = "#377EB8", Adenoma = "#E6AB02",
      CRC = "#E41A1C", `colorectal carcinoma` = "#E41A1C"), drop = FALSE) +
    labs(title = expression(italic("F. nucleatum") * " baseline detection varies by profiler"),
         subtitle = "Unspiked matched baselines; labels are detected samples / samples",
         x = NULL, y = "Baseline prevalence", fill = "Condition") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  figure_sources$fnuc_baseline_detection <- fnuc_prevalence
}

if (!is.null(fnuc_classes)) {
  require_columns(fnuc_classes, c("cohort", "analysis_population", "condition", "profiler",
    "nominal_target_fraction", "recovery_class", "samples", "fraction"),
    "fnuc_recovery_classes.tsv")
  fnuc_classes$nominal_target_fraction <- numeric_column(
    fnuc_classes, "nominal_target_fraction", "fnuc_recovery_classes.tsv", 0)
  fnuc_classes$fraction <- numeric_column(
    fnuc_classes, "fraction", "fnuc_recovery_classes.tsv", 0, 1)
  fnuc_classes <- decorate(fnuc_classes)
  fnuc_classes$Dose <- factor(
    percent(fnuc_classes$nominal_target_fraction, accuracy = .001),
    levels = percent(sort(unique(fnuc_classes$nominal_target_fraction)), accuracy = .001))
  fnuc_classes$recovery_class <- factor(fnuc_classes$recovery_class,
    levels = c("Poor / missed", "Average", "Good"))
  fnuc_classes <- fnuc_classes[fnuc_classes$Population == "Community" &
    fnuc_classes$nominal_target_fraction %in% c(0.0001, 0.0005, 0.001), , drop = FALSE]
  fnuc_classes$Dose <- factor(fnuc_classes$Dose,
    levels = percent(c(0.0001, 0.0005, 0.001), accuracy = .001))
  fnuc_classes$Condition <- factor(fnuc_classes$Condition,
    levels = c("Control", "Adenoma", "CRC", "colorectal carcinoma"))
  # Pool cohorts for the compact poster-style comparison.  Fractions are
  # recomputed from pooled sample counts, not averaged across cohorts.
  fnuc_classes <- group_apply(fnuc_classes,
    c("analysis_population", "condition", "profiler", "nominal_target_fraction", "recovery_class"),
    function(x) data.frame(samples = sum(x$samples)))
  totals <- group_apply(fnuc_classes,
    c("analysis_population", "condition", "profiler", "nominal_target_fraction"),
    function(x) data.frame(total_samples = sum(x$samples)))
  fnuc_classes <- merge(fnuc_classes, totals,
    by = c("analysis_population", "condition", "profiler", "nominal_target_fraction"),
    all.x = TRUE, sort = FALSE)
  fnuc_classes$fraction <- fnuc_classes$samples / fnuc_classes$total_samples
  fnuc_classes <- decorate(fnuc_classes)
  fnuc_classes$Dose <- factor(
    percent(fnuc_classes$nominal_target_fraction, accuracy = .001),
    levels = percent(c(0.0001, 0.0005, 0.001), accuracy = .001))
  plots$fnuc_recovery_classes <- ggplot(
    fnuc_classes, aes(Condition, fraction, fill = recovery_class)
  ) +
    geom_col(width = .78) +
    geom_text(aes(label = ifelse(fraction >= .08, paste0(round(100*fraction), "%"), "")),
              position = position_stack(vjust = .5), size = 2.7) +
    facet_grid(Profiler ~ Dose) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = c("Poor / missed" = "#D55E00", Average = "#E69F00",
                                 Good = "#009E73"), drop = FALSE) +
    labs(title = expression(italic("F. nucleatum") * " recovery across the spike series"),
         subtitle = "Good: <=10% relative error; Average: 10-50%; Poor/missed: >50% or undetected",
         x = "Implanted target fraction", y = "Samples", fill = "Recovery") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  figure_sources$fnuc_recovery_classes <- fnuc_classes
}

if (!is.null(biomarker_thresholds)) {
  require_columns(biomarker_thresholds, c("cohort", "analysis_population",
    "target_label", "profiler", "minimum_called_fraction", "ever_called", "tested_doses"),
    "biomarker_detection_thresholds.tsv")
  biomarker_thresholds$minimum_called_fraction <- suppressWarnings(
    as.numeric(biomarker_thresholds$minimum_called_fraction))
  biomarker_thresholds$ever_called <- as.integer(biomarker_thresholds$ever_called)
  biomarker_thresholds <- decorate(biomarker_thresholds)
  biomarker_thresholds$threshold_label <- ifelse(
    biomarker_thresholds$ever_called == 1,
    percent(biomarker_thresholds$minimum_called_fraction, accuracy = .001), "NR")
  biomarker_thresholds$threshold_log10 <- ifelse(
    biomarker_thresholds$ever_called == 1,
    log10(biomarker_thresholds$minimum_called_fraction), NA_real_)
  biomarker_thresholds$label_colour <- ifelse(
    is.na(biomarker_thresholds$threshold_log10) |
      biomarker_thresholds$threshold_log10 > -3, "dark", "light")
  plots$paired_biomarker_detection_threshold <- ggplot(
    biomarker_thresholds, aes(Cohort, target_label, fill = threshold_log10)
  ) +
    geom_tile(colour = "white", linewidth = .35) +
    geom_text(aes(label = threshold_label, colour = label_colour), size = 2.8) +
    facet_grid(Population ~ Profiler, scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colours = c("#00204D", "#414D6B", "#6C6E72", "#9C9165", "#C8B35B", "#F1D54A"),
      na.value = "#BDBDBD",
      labels = function(x) percent(10^x, accuracy = .001)) +
    scale_colour_manual(values = c(dark = "black", light = "white"), guide = "none") +
    labs(title = "Minimum perturbation required for paired biomarker discovery",
         subtitle = "Spiked versus matched baseline, pooled across conditions; BH q <= 0.05; NR = not reached",
         x = NULL, y = "Implanted target", fill = "Minimum fraction")
  figure_sources$paired_biomarker_detection_threshold <- biomarker_thresholds
}

dimensions <- list(
  response_operator_crosstalk = c(12, 8),
  response_operator_performance_map = c(11, 7),
  superposition_error_comparison = c(11, 9),
  superposition_improvement_by_dose = c(11, 8),
  measurement_reliability_landscape = c(12, 8),
  measurement_reliability_transport = c(11, 7),
  detection_instability_by_dose = c(14, 8),
  quantitative_error_by_dose = c(14, 8),
  sentinel_panel_saturation = c(12, 8),
  heldout_reliability_prediction = c(12, 8),
  heldout_operator_generalization = c(11, 7),
  cohort_holdout_reliability_summary = c(10, 8),
  biomarker_replication_by_reliability = c(10, 7),
  implanted_target_recovery_map = c(15, 10),
  fnuc_baseline_detection = c(12, 5.5),
  fnuc_recovery_classes = c(14, 11),
  paired_biomarker_detection_threshold = c(12, 10)
)
for (name in names(plots)) {
  size <- dimensions[[name]]
  ggsave(file.path(outdir, "figures", paste0(name, ".pdf")), plots[[name]],
         width = size[[1L]], height = size[[2L]])
  ggsave(file.path(outdir, "figures", paste0(name, ".png")), plots[[name]],
         width = size[[1L]], height = size[[2L]], dpi = 300, bg = "white")
  write_tsv(figure_sources[[name]],
            file.path(outdir, "figure_source", paste0(name, ".tsv")))
}
if ("measurement_reliability_summary" %in% names(figure_sources)) {
  write_tsv(figure_sources$measurement_reliability_summary,
            file.path(outdir, "tables", "measurement_reliability_summary.tsv"))
}
if ("measurement_reliability_estimability" %in% names(figure_sources)) {
  write_tsv(figure_sources$measurement_reliability_estimability,
            file.path(outdir, "tables", "measurement_reliability_estimability.tsv"))
}

for (d in input_tables) {
  source <- attr(d, "source_path")
  file.copy(source, file.path(outdir, "tables", basename(source)), overwrite = TRUE)
}

available <- data.frame(
  analysis = c("response_operator", "superposition", "reliability_certificates",
               "detection_response", "quantitative_response", "heldout_validation",
               "cohort_holdout_summary", "panel_saturation", "biomarker_replication"),
  status = c("GENERATED", "GENERATED", "GENERATED",
             if (is.null(detection)) "NOT_PROVIDED" else "GENERATED",
             if (is.null(quantitative)) "NOT_PROVIDED" else "GENERATED",
             if (is.null(heldout)) "NOT_PROVIDED" else "GENERATED",
             if (is.null(cohort_holdout)) "NOT_PROVIDED" else "GENERATED",
             if (is.null(panel)) "NOT_PROVIDED" else "GENERATED",
             if (is.null(biomarker)) "NOT_PROVIDED" else "GENERATED"),
  stringsAsFactors = FALSE
)
write_tsv(available, file.path(outdir, "diagnostics", "analysis_availability.tsv"))
diagnostics <- data.frame(
  metric = c("report_status", "figures", "operator_rows", "superposition_rows",
             "certificate_rows", "input_tables"),
  value = c(report_status, length(plots), nrow(operator), nrow(superposition),
            nrow(certificates), length(input_tables)), stringsAsFactors = FALSE
)
write_tsv(diagnostics, file.path(outdir, "diagnostics", "report_diagnostics.tsv"))

guide <- c(
  "# Target-agnostic perturbation-response report", "",
  "The response-operator panels separate correct target gain from off-diagonal bystander cross-talk.",
  "The superposition panels test whether responses learned from independent single-taxon spikes",
  "predict community mixtures better than the exact composition-only counterfactual.",
  "Reliability certificates are learned while each reported feature is a bystander; they do not",
  "require knowing which taxon may perturb a future sample.",
  "Held-out and panel-saturation figures are generated only when their source tables are supplied.",
  "External biomarker replication is an optional validation endpoint and never enters the score.", "",
  paste0("Report status: ", report_status, ".")
)
writeLines(guide, file.path(outdir, "FIGURE_GUIDE.md"))
if (report_status == "DEVELOPMENT_ONLY") {
  writeLines(c("status=DEVELOPMENT_ONLY", "use_for_manuscript=NO",
               "thresholds_frozen=NO"), file.path(outdir, "DEVELOPMENT_ONLY.txt"))
}
manifest <- data.frame(
  field = c("report", "status", "input_root", "created_at"),
  value = c("target_agnostic_perturbation_response", report_status,
            normalizePath(input_root), format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(outdir, "provenance", "report_manifest.tsv"))

checksum <- file.path(outdir, "provenance", "perturbation_response_report.sha256")
outputs <- list.files(outdir, recursive = TRUE, full.names = TRUE)
files <- c(input_paths, outputs[!grepl("SUCCESS$|perturbation_response_report\\.sha256$", outputs)])
rc <- system2("sha256sum", files, stdout = checksum)
if (!identical(rc, 0L)) stop("Could not seal perturbation-response report")
writeLines(c("report\ttarget_agnostic_perturbation_response",
             paste0("report_status\t", report_status),
             paste0("figures\t", length(plots)), "status\tPASS"),
           file.path(outdir, "SUCCESS"))
message("[PASS] Target-agnostic perturbation-response report: ", normalizePath(outdir))
