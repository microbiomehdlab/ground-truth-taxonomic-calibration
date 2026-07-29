#!/usr/bin/env Rscript

# Spike-informed calibratability analysis using linear models and GAMs.
# Primary goals:
#   1) Visualise expected-versus-observed abundance response curves.
#   2) Compare raw identity, linear, and GAM inverse calibration models
#      using grouped cross-validation by base sample.
#   3) Summarise profiler–taxon combinations for which held-out prediction
#      improves over the uncorrected abundance scale.
#
# Required packages: mgcv, ggplot2
# Optional package: patchwork (for a composite figure)

parse_args <- function(args) {
  out <- list(
    input = "RUNS/spike_metrics/community/target_member_errors_with_condition.csv",
    outdir = "RUNS/manuscript_figures/spike_calibratability_gam",
    folds = 5L,
    repeats = 5L,
    seed = 1L,
    k = 5L,
    min_rows = 80L,
    min_detected = 40L,
    min_base_samples = 8L,
    min_detection_rate = 0.50,
    strong_improvement = 0.10,
    strong_max_log_error = 0.25,
    pseudocount = NA_real_,
    focus_targets = "Bacteroides fragilis,Fusobacterium nucleatum subsp. nucleatum,Parvimonas micra,Dialister pneumosintes",
    transfer_tests = TRUE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(args)) stop("Missing value for --", key)
    val <- args[[i + 1L]]
    key2 <- gsub("-", "_", key)
    if (!key2 %in% names(out)) stop("Unknown option --", key)
    if (key2 %in% c(
      "folds", "repeats", "seed", "k", "min_rows", "min_detected",
      "min_base_samples"
    )) {
      val <- as.integer(val)
    } else if (key2 %in% c(
      "pseudocount", "min_detection_rate", "strong_improvement",
      "strong_max_log_error"
    )) {
      val <- as.numeric(val)
    } else if (key2 == "transfer_tests") {
      val <- tolower(val) %in% c("true", "t", "1", "yes", "y")
    }
    out[[key2]] <- val
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (args$folds < 2L) stop("--folds must be at least 2.")
if (args$repeats < 1L) stop("--repeats must be at least 1.")
if (args$k < 3L) stop("--k must be at least 3.")
if (args$min_rows < 1L) stop("--min-rows must be at least 1.")
if (args$min_detected < 1L) stop("--min-detected must be at least 1.")
if (args$min_base_samples < 2L) stop("--min-base-samples must be at least 2.")
if (!is.finite(args$min_detection_rate) ||
    args$min_detection_rate < 0 ||
    args$min_detection_rate > 1) {
  stop("--min-detection-rate must lie between 0 and 1.")
}
if (!is.finite(args$strong_improvement) || args$strong_improvement < 0) {
  stop("--strong-improvement must be non-negative.")
}
if (!is.finite(args$strong_max_log_error) ||
    args$strong_max_log_error <= 0) {
  stop("--strong-max-log-error must be positive.")
}

if (!file.exists(args$input)) stop("Input file not found: ", args$input)
if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
message("[INFO] Reading: ", args$input)
d <- read.csv(args$input, check.names = FALSE, stringsAsFactors = FALSE)

required <- c("base_id", "tool", "member_taxon", "observed", "expected", "baseline",
              "detected", "Target_Condition", "Study", "spike_fraction_total")
missing_cols <- setdiff(required, names(d))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# Standardise types and retain canonical target rows.
d$observed <- suppressWarnings(as.numeric(d$observed))
d$expected <- suppressWarnings(as.numeric(d$expected))
d$baseline <- suppressWarnings(as.numeric(d$baseline))
d$spike_fraction_total <- suppressWarnings(as.numeric(d$spike_fraction_total))
d$detected <- as.logical(d$detected)
d <- d[is.finite(d$observed) & is.finite(d$expected) & d$expected >= 0 & d$observed >= 0, , drop = FALSE]
d <- d[!is.na(d$member_taxon) & nzchar(d$member_taxon), , drop = FALSE]

# Derive a reproducible pseudocount from the data unless supplied.
nonzero <- c(d$observed[d$observed > 0], d$expected[d$expected > 0])
if (!length(nonzero)) stop("No positive observed or expected abundances found.")
pc <- args$pseudocount
if (!is.finite(pc) || pc <= 0) pc <- min(nonzero, na.rm = TRUE) / 2
message("[INFO] Pseudocount: ", format(pc, scientific = TRUE))

d$log_expected <- log10(d$expected + pc)
d$log_observed <- log10(d$observed + pc)
d$pair_id <- interaction(d$tool, d$member_taxon, drop = TRUE, sep = " || ")
d$tool_label <- ifelse(d$tool == "kraken2_bracken", "Kraken2 + Bracken",
                       ifelse(grepl("metaphlan", d$tool, ignore.case = TRUE), "MetaPhlAn 4", d$tool))

# Paper label codes used throughout the manuscript/table.
# The canonical taxon names are retained in all output tables; these short labels
# are used only for figure display.
taxon_label_map <- c(
  "Bacteroides fragilis" = "Bfrag",
  "Clostridium symbiosum" = "Csym",
  "Dialister pneumosintes" = "Dpne",
  "Fusobacterium nucleatum subsp. nucleatum" = "Fnuc",
  "Hungatella hathewayi" = "Hhat",
  "Parvimonas micra" = "Pmic",
  "Peptostreptococcus anaerobius" = "Pana",
  "Peptostreptococcus stomatis" = "Psto",
  "Porphyromonas asaccharolytica" = "Porp",
  "Prevotella intermedia" = "Pint"
)

# Shared manuscript order used in Figures 3–5 and the supplementary figures.
taxon_order <- c(
  "Bacteroides fragilis",
  "Clostridium symbiosum",
  "Dialister pneumosintes",
  "Fusobacterium nucleatum subsp. nucleatum",
  "Hungatella hathewayi",
  "Parvimonas micra",
  "Peptostreptococcus anaerobius",
  "Peptostreptococcus stomatis",
  "Porphyromonas asaccharolytica",
  "Prevotella intermedia"
)
d$taxon_label <- unname(taxon_label_map[d$member_taxon])
d$taxon_label[is.na(d$taxon_label)] <- d$member_taxon[is.na(d$taxon_label)]

# Detection summaries are kept separate from abundance-response modelling.
detection_summary <- aggregate(
  as.numeric(d$detected),
  by = list(tool = d$tool_label, taxon = d$member_taxon, fraction = d$spike_fraction_total),
  FUN = function(x) c(n = length(x), detected = sum(x, na.rm = TRUE), rate = mean(x, na.rm = TRUE))
)
detection_summary <- do.call(data.frame, detection_summary)
names(detection_summary)[4:6] <- c("n", "n_detected", "detection_rate")
write.csv(detection_summary, file.path(args$outdir, "detection_summary.csv"), row.names = FALSE)

# Only positive detections are used for quantitative abundance-response models.
positive_flag <- !is.na(d$detected) & d$detected &
  d$observed > 0 & d$expected > 0
q <- d[positive_flag, , drop = FALSE]

# Construct a complete eligibility table before modelling so that profiler--taxon
# combinations with insufficient detection or modelling data remain visible and
# are classified as not assessable rather than silently disappearing.
effective_min_base_samples <- max(args$min_base_samples, args$folds)
all_pairs <- split(d, d$pair_id)

summarise_pair_eligibility <- function(z) {
  detected_flag <- !is.na(z$detected) & z$detected
  positive <- detected_flag & z$observed > 0 & z$expected > 0
  positive_base <- z$base_id[positive & !is.na(z$base_id)]

  detection_rate <- mean(z$detected, na.rm = TRUE)
  if (!is.finite(detection_rate)) detection_rate <- NA_real_

  data.frame(
    pair_id = as.character(z$pair_id[1]),
    tool = z$tool_label[1],
    taxon = z$member_taxon[1],
    n_all = nrow(z),
    n_detected = sum(detected_flag),
    detection_rate_all = detection_rate,
    n_positive = sum(positive),
    n_base_positive = length(unique(positive_base)),
    stringsAsFactors = FALSE
  )
}

pair_eligibility <- do.call(
  rbind,
  lapply(all_pairs, summarise_pair_eligibility)
)
rownames(pair_eligibility) <- NULL

pair_eligibility$input_assessable <- with(
  pair_eligibility,
  is.finite(detection_rate_all) &
    detection_rate_all >= args$min_detection_rate &
    n_detected >= args$min_detected &
    n_positive >= args$min_rows &
    n_base_positive >= effective_min_base_samples
)

make_eligibility_reason <- function(i) {
  reasons <- character()

  if (!is.finite(pair_eligibility$detection_rate_all[i])) {
    reasons <- c(reasons, "Detection rate unavailable")
  } else if (pair_eligibility$detection_rate_all[i] <
             args$min_detection_rate) {
    reasons <- c(
      reasons,
      sprintf(
        "Detection rate below %.2f",
        args$min_detection_rate
      )
    )
  }

  if (pair_eligibility$n_detected[i] < args$min_detected) {
    reasons <- c(
      reasons,
      sprintf(
        "Fewer than %d detected observations",
        args$min_detected
      )
    )
  }

  if (pair_eligibility$n_positive[i] < args$min_rows) {
    reasons <- c(
      reasons,
      sprintf(
        "Fewer than %d positive-abundance observations",
        args$min_rows
      )
    )
  }

  if (pair_eligibility$n_base_positive[i] <
      effective_min_base_samples) {
    reasons <- c(
      reasons,
      sprintf(
        "Fewer than %d biological samples with positive observations",
        effective_min_base_samples
      )
    )
  }

  if (!length(reasons)) "Eligible for grouped cross-validation" else
    paste(reasons, collapse = "; ")
}

pair_eligibility$eligibility_reason <- vapply(
  seq_len(nrow(pair_eligibility)),
  make_eligibility_reason,
  character(1)
)
pair_eligibility$cv_success <- FALSE

eligible_pair_ids <- pair_eligibility$pair_id[
  pair_eligibility$input_assessable
]
eligible_q <- q[
  as.character(q$pair_id) %in% eligible_pair_ids,
  ,
  drop = FALSE
]
eligible_q$pair_id <- droplevels(eligible_q$pair_id)
pairs <- split(
  eligible_q,
  eligible_q$pair_id,
  drop = TRUE
)

threshold_table <- data.frame(
  parameter = c(
    "min_detection_rate",
    "min_detected",
    "min_positive_observations",
    "min_base_samples",
    "strong_min_log_error_improvement",
    "strong_max_median_abs_log_error",
    "grouped_cv_folds",
    "grouped_cv_repeats"
  ),
  value = c(
    args$min_detection_rate,
    args$min_detected,
    args$min_rows,
    effective_min_base_samples,
    args$strong_improvement,
    args$strong_max_log_error,
    args$folds,
    args$repeats
  ),
  interpretation = c(
    "Minimum overall detection rate required for assessment",
    "Minimum number of detected observations required for assessment",
    "Minimum number of positive-abundance observations used for modelling",
    "Minimum number of original biological samples contributing positive observations",
    "Strictly greater improvement required for strong evidence",
    "Strictly lower residual median absolute log10 error required for strong evidence",
    "Number of grouped cross-validation folds",
    "Number of repeated grouped cross-validation runs"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  threshold_table,
  file.path(args$outdir, "calibratability_thresholds.csv"),
  row.names = FALSE
)

make_group_folds <- function(base_id, k, seed) {
  ids <- unique(base_id)
  set.seed(seed)
  ids <- sample(ids)
  fold_id <- rep(seq_len(min(k, length(ids))), length.out = length(ids))
  names(fold_id) <- ids
  unname(fold_id[base_id])
}

metric_row <- function(truth_log, pred_log, truth, pred, model, repeat_id, fold_id, tool, taxon) {
  ok <- is.finite(truth_log) & is.finite(pred_log) & is.finite(truth) & is.finite(pred)
  truth_log <- truth_log[ok]; pred_log <- pred_log[ok]; truth <- truth[ok]; pred <- pmax(pred[ok], 0)
  rel <- abs(pred - truth) / pmax(truth, pc)
  data.frame(
    tool = tool,
    taxon = taxon,
    model = model,
    cv_repeat = repeat_id,
    fold = fold_id,
    n_test = length(truth),
    median_abs_log_error = median(abs(pred_log - truth_log), na.rm = TRUE),
    rmse_log = sqrt(mean((pred_log - truth_log)^2, na.rm = TRUE)),
    median_abs_relative_error = median(rel, na.rm = TRUE),
    within_10pct = mean(rel <= 0.10, na.rm = TRUE),
    within_50pct = mean(rel <= 0.50, na.rm = TRUE),
    spearman = suppressWarnings(cor(truth, pred, method = "spearman", use = "complete.obs")),
    stringsAsFactors = FALSE
  )
}

cv_predictions <- list()
cv_metrics <- list()
fit_diagnostics <- list()
idx <- 1L; midx <- 1L; didx <- 1L

message(
  "[INFO] Evaluating ",
  length(pairs),
  " eligible profiler–taxon combinations (",
  nrow(pair_eligibility),
  " total combinations)"
)

for (nm in names(pairs)) {
  z <- pairs[[nm]]
  tool <- z$tool_label[1]
  taxon <- z$member_taxon[1]
  pair_metric_start <- midx

  # Forward model for curve diagnostics: observed as a function of expected.
  k_use <- min(args$k, max(3L, length(unique(z$log_expected)) - 1L))
  forward_gam <- try(mgcv::gam(log_observed ~ s(log_expected, k = k_use, bs = "cr"),
                               data = z, method = "REML"), silent = TRUE)
  if (!inherits(forward_gam, "try-error")) {
    sm <- summary(forward_gam)
    edf <- if (nrow(sm$s.table)) sm$s.table[1, "edf"] else NA_real_
    fit_diagnostics[[didx]] <- data.frame(
      tool = tool, taxon = taxon, n = nrow(z), n_base = length(unique(z$base_id)),
      detection_rate_all = mean(d[d$pair_id == nm, "detected"], na.rm = TRUE),
      forward_gam_edf = edf,
      forward_deviance_explained = sm$dev.expl,
      stringsAsFactors = FALSE
    )
    didx <- didx + 1L
  }

  for (r in seq_len(args$repeats)) {
    folds <- make_group_folds(z$base_id, args$folds, args$seed + 1000L * r + didx)
    for (f in sort(unique(folds))) {
      train <- z[folds != f, , drop = FALSE]
      test <- z[folds == f, , drop = FALSE]
      if (nrow(train) < 20 || nrow(test) < 5) next

      # Inverse calibration: estimate expected abundance from observed abundance.
      lm_fit <- try(lm(log_expected ~ log_observed, data = train), silent = TRUE)
      gam_fit <- try(mgcv::gam(log_expected ~ s(log_observed, k = k_use, bs = "cr"),
                               data = train, method = "REML"), silent = TRUE)

      pred_identity_log <- test$log_observed
      pred_lm_log <- if (!inherits(lm_fit, "try-error")) as.numeric(predict(lm_fit, newdata = test)) else rep(NA_real_, nrow(test))
      pred_gam_log <- if (!inherits(gam_fit, "try-error")) as.numeric(predict(gam_fit, newdata = test)) else rep(NA_real_, nrow(test))

      pred_identity <- pmax(10^pred_identity_log - pc, 0)
      pred_lm <- pmax(10^pred_lm_log - pc, 0)
      pred_gam <- pmax(10^pred_gam_log - pc, 0)

      pred_block <- rbind(
        data.frame(test[c("base_id", "sample_id", "Study", "Target_Condition", "spike_fraction_total")],
                   tool = tool, taxon = taxon, cv_repeat = r, fold = f, model = "Raw identity",
                   expected = test$expected, observed = test$observed, predicted_expected = pred_identity),
        data.frame(test[c("base_id", "sample_id", "Study", "Target_Condition", "spike_fraction_total")],
                   tool = tool, taxon = taxon, cv_repeat = r, fold = f, model = "Linear",
                   expected = test$expected, observed = test$observed, predicted_expected = pred_lm),
        data.frame(test[c("base_id", "sample_id", "Study", "Target_Condition", "spike_fraction_total")],
                   tool = tool, taxon = taxon, cv_repeat = r, fold = f, model = "GAM",
                   expected = test$expected, observed = test$observed, predicted_expected = pred_gam)
      )
      cv_predictions[[idx]] <- pred_block
      idx <- idx + 1L

      cv_metrics[[midx]] <- rbind(
        metric_row(test$log_expected, pred_identity_log, test$expected, pred_identity, "Raw identity", r, f, tool, taxon),
        metric_row(test$log_expected, pred_lm_log, test$expected, pred_lm, "Linear", r, f, tool, taxon),
        metric_row(test$log_expected, pred_gam_log, test$expected, pred_gam, "GAM", r, f, tool, taxon)
      )
      midx <- midx + 1L
    }
  }

  eligibility_index <- match(
    as.character(nm),
    pair_eligibility$pair_id
  )
  if (!is.na(eligibility_index)) {
    pair_eligibility$cv_success[eligibility_index] <-
      midx > pair_metric_start
  }
}

pair_eligibility$assessment_reason <- pair_eligibility$eligibility_reason
cv_failed <- pair_eligibility$input_assessable &
  !pair_eligibility$cv_success
pair_eligibility$assessment_reason[cv_failed] <-
  "No valid grouped cross-validation estimate was obtained"

write.csv(
  pair_eligibility,
  file.path(args$outdir, "calibratability_input_eligibility.csv"),
  row.names = FALSE
)

if (!length(cv_metrics)) stop("No profiler–taxon combination met the modelling requirements.")
cv_metrics <- do.call(rbind, cv_metrics)
cv_predictions <- do.call(rbind, cv_predictions)
fit_diagnostics <- if (length(fit_diagnostics)) do.call(rbind, fit_diagnostics) else data.frame()
write.csv(cv_metrics, file.path(args$outdir, "grouped_cv_metrics_by_fold.csv"), row.names = FALSE)
write.csv(cv_predictions, file.path(args$outdir, "grouped_cv_predictions.csv"), row.names = FALSE)
if (nrow(fit_diagnostics)) write.csv(fit_diagnostics, file.path(args$outdir, "forward_gam_diagnostics.csv"), row.names = FALSE)

# Aggregate cross-validation performance.
metric_cols <- c("median_abs_log_error", "rmse_log", "median_abs_relative_error", "within_10pct", "within_50pct", "spearman")
agg <- aggregate(cv_metrics[metric_cols], by = cv_metrics[c("tool", "taxon", "model")], FUN = median, na.rm = TRUE)
raw <- agg[agg$model == "Raw identity", c("tool", "taxon", "median_abs_log_error", "median_abs_relative_error")]
names(raw)[3:4] <- c("raw_median_abs_log_error", "raw_median_abs_relative_error")
agg <- merge(agg, raw, by = c("tool", "taxon"), all.x = TRUE)
agg$log_error_improvement <- agg$raw_median_abs_log_error - agg$median_abs_log_error
agg$relative_error_improvement <- agg$raw_median_abs_relative_error - agg$median_abs_relative_error
write.csv(agg, file.path(args$outdir, "grouped_cv_model_comparison_summary.csv"), row.names = FALSE)

# Calibratability summary based on held-out GAM improvement, residual error,
# detection frequency, and the availability of sufficient modelling data.
gam_metrics <- agg[agg$model == "GAM", , drop = FALSE]
gam_sum <- merge(
  pair_eligibility,
  gam_metrics,
  by = c("tool", "taxon"),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(fit_diagnostics)) {
  diagnostic_columns <- c(
    "tool", "taxon",
    "forward_gam_edf",
    "forward_deviance_explained"
  )
  gam_sum <- merge(
    gam_sum,
    fit_diagnostics[diagnostic_columns],
    by = c("tool", "taxon"),
    all.x = TRUE,
    sort = FALSE
  )
}

valid_gam_metrics <- with(
  gam_sum,
  input_assessable &
    cv_success &
    is.finite(log_error_improvement) &
    is.finite(median_abs_log_error)
)

gam_sum$calibratability_class <- "Not assessable"
gam_sum$calibratability_code <- "ND"
gam_sum$classification_reason <- gam_sum$assessment_reason

missing_finite_gam <- gam_sum$input_assessable &
  gam_sum$cv_success &
  !valid_gam_metrics
gam_sum$classification_reason[missing_finite_gam] <-
  "No finite held-out GAM performance estimate was obtained"

strong_index <- valid_gam_metrics &
  gam_sum$log_error_improvement > args$strong_improvement &
  gam_sum$median_abs_log_error < args$strong_max_log_error

conditional_index <- valid_gam_metrics &
  !strong_index &
  gam_sum$log_error_improvement > 0

little_index <- valid_gam_metrics &
  !strong_index &
  !conditional_index

gam_sum$calibratability_class[strong_index] <- "Strong evidence"
gam_sum$calibratability_code[strong_index] <- "S"
gam_sum$classification_reason[strong_index] <- sprintf(
  "Improvement > %.2f and residual median absolute log10 error < %.2f",
  args$strong_improvement,
  args$strong_max_log_error
)

gam_sum$calibratability_class[conditional_index] <-
  "Conditional evidence"
gam_sum$calibratability_code[conditional_index] <- "C"
gam_sum$classification_reason[conditional_index] <-
  "Positive held-out improvement, but one or both strong-evidence criteria were not met"

gam_sum$calibratability_class[little_index] <- "Little evidence"
gam_sum$calibratability_code[little_index] <- "L"
gam_sum$classification_reason[little_index] <-
  "No positive held-out improvement over the raw identity prediction"

gam_sum$min_detection_rate_threshold <- args$min_detection_rate
gam_sum$min_detected_threshold <- args$min_detected
gam_sum$min_positive_threshold <- args$min_rows
gam_sum$min_base_samples_threshold <- effective_min_base_samples
gam_sum$strong_improvement_threshold <- args$strong_improvement
gam_sum$strong_max_log_error_threshold <- args$strong_max_log_error

gam_sum <- gam_sum[
  order(gam_sum$tool, gam_sum$taxon),
  ,
  drop = FALSE
]

write.csv(
  gam_sum,
  file.path(args$outdir, "calibratability_summary.csv"),
  row.names = FALSE
)

# Optional cross-study transfer test.
transfer_rows <- list(); tidx <- 1L
if (isTRUE(args$transfer_tests) && length(unique(q$Study)) >= 2) {
  for (nm in names(pairs)) {
    z <- pairs[[nm]]
    if (nrow(z) < args$min_rows) next
    studies <- unique(z$Study)
    if (length(studies) < 2) next
    k_use <- min(args$k, max(3L, length(unique(z$log_observed)) - 1L))
    for (train_study in studies) {
      test_studies <- setdiff(studies, train_study)
      train <- z[z$Study == train_study, , drop = FALSE]
      if (nrow(train) < 30) next
      fit <- try(mgcv::gam(log_expected ~ s(log_observed, k = k_use, bs = "cr"), data = train, method = "REML"), silent = TRUE)
      if (inherits(fit, "try-error")) next
      for (test_study in test_studies) {
        test <- z[z$Study == test_study, , drop = FALSE]
        if (nrow(test) < 10) next
        pred_log <- as.numeric(predict(fit, newdata = test))
        pred <- pmax(10^pred_log - pc, 0)
        raw_log <- test$log_observed
        raw_pred <- test$observed
        transfer_rows[[tidx]] <- rbind(
          transform(metric_row(test$log_expected, raw_log, test$expected, raw_pred, "Raw identity", 1, 1,
                               z$tool_label[1], z$member_taxon[1]), train_study = train_study, test_study = test_study),
          transform(metric_row(test$log_expected, pred_log, test$expected, pred, "GAM", 1, 1,
                               z$tool_label[1], z$member_taxon[1]), train_study = train_study, test_study = test_study)
        )
        tidx <- tidx + 1L
      }
    }
  }
}
transfer <- if (length(transfer_rows)) do.call(rbind, transfer_rows) else data.frame()
if (nrow(transfer)) write.csv(transfer, file.path(args$outdir, "cross_study_transfer_metrics.csv"), row.names = FALSE)

# ----------------------------- Figures ---------------------------------
library(ggplot2)
theme_manuscript <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom", plot.title.position = "plot")

percent_axis_label <- function(x) {
  p <- as.numeric(x) * 100
  ifelse(
    p < 0.01,
    sprintf("%.3f%%", p),
    ifelse(p < 1, sprintf("%.2f%%", p), paste0(formatC(p, format = "fg", digits = 3, drop0trailing = TRUE), "%"))
  )
}

focus <- trimws(strsplit(args$focus_targets, ",", fixed = TRUE)[[1]])
focus <- focus[focus %in% unique(q$member_taxon)]
if (!length(focus)) focus <- head(unique(q$member_taxon), 4)
plot_dat <- q[q$member_taxon %in% focus, , drop = FALSE]
focus_order <- taxon_order[taxon_order %in% focus]
plot_dat$taxon_label <- factor(
  plot_dat$taxon_label,
  levels = unname(taxon_label_map[focus_order])
)

# Panel A: sample-level expected vs observed abundance.
# The GAM is fitted on log10-transformed abundance, matching the model used in the
# calibratability analysis, and predictions are transformed back to the original
# abundance scale for plotting.
build_forward_curve <- function(z, k_max = 5L, n_grid = 200L) {
  ux <- length(unique(z$log_expected[is.finite(z$log_expected)]))
  k_use <- min(k_max, max(3L, ux - 1L))
  fit <- try(
    mgcv::gam(
      log_observed ~ s(log_expected, k = k_use, bs = "cr"),
      data = z,
      method = "REML"
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)

  grid <- data.frame(
    log_expected = seq(min(z$log_expected, na.rm = TRUE),
                       max(z$log_expected, na.rm = TRUE),
                       length.out = n_grid)
  )
  pr <- predict(fit, newdata = grid, se.fit = TRUE)
  crit <- qnorm(0.975)

  data.frame(
    tool_label = z$tool_label[1],
    taxon_label = z$taxon_label[1],
    expected = pmax(10^grid$log_expected - pc, pc / 10),
    fitted = pmax(10^as.numeric(pr$fit) - pc, pc / 10),
    lower = pmax(10^(as.numeric(pr$fit) - crit * as.numeric(pr$se.fit)) - pc, pc / 10),
    upper = pmax(10^(as.numeric(pr$fit) + crit * as.numeric(pr$se.fit)) - pc, pc / 10),
    stringsAsFactors = FALSE
  )
}

curve_split <- split(plot_dat, interaction(plot_dat$tool_label, plot_dat$taxon_label, drop = TRUE))
forward_curves <- lapply(curve_split, build_forward_curve, k_max = args$k)
forward_curves <- forward_curves[!vapply(forward_curves, is.null, logical(1))]
forward_curves <- if (length(forward_curves)) do.call(rbind, forward_curves) else data.frame()

pA <- ggplot(plot_dat, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4) +
  geom_point(alpha = 0.22, size = 0.7) +
  {if (nrow(forward_curves)) geom_ribbon(
    data = forward_curves,
    aes(x = expected, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    alpha = 0.12
  )} +
  {if (nrow(forward_curves)) geom_line(
    data = forward_curves,
    aes(x = expected, y = fitted),
    inherit.aes = FALSE,
    linewidth = 0.7
  )} +
  scale_x_log10(
    labels = percent_axis_label
  ) +
  scale_y_log10(
    labels = percent_axis_label
  ) +
  facet_grid(tool_label ~ taxon_label, scales = "free") +
  labs(title = "A. Expected and observed abundance for representative taxa",
       subtitle = "Points are samples; dashed line is identity; solid line is the GAM fit and ribbon is its 95% confidence interval",
       x = "Expected total post-spike abundance (%)",
       y = "Observed total post-spike abundance (%)") + theme_manuscript

pA_individual <- pA + labs(title = NULL, subtitle = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
ggsave(file.path(args$outdir, "panel_A_expected_vs_observed_gam.pdf"), pA_individual, width = 13, height = 6.8, device = cairo_pdf)
ggsave(file.path(args$outdir, "panel_A_expected_vs_observed_gam.png"), pA_individual, width = 13, height = 6.8, dpi = 450)

# Panel B: model comparison.
model_order <- c("Raw identity", "Linear", "GAM")
agg$model <- factor(agg$model, levels = model_order)
pB <- ggplot(agg, aes(x = model, y = median_abs_log_error, group = interaction(tool, taxon))) +
  geom_line(alpha = 0.25) + geom_point(size = 1.5) +
  facet_wrap(~ tool, scales = "free_y") +
  labs(title = "B. Grouped cross-validated model performance",
       subtitle = "Lower held-out median absolute log error is better",
       x = NULL, y = "Median absolute log10 error") + theme_manuscript

pB_individual <- pB + labs(title = NULL, subtitle = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
ggsave(file.path(args$outdir, "panel_B_model_comparison.pdf"), pB_individual, width = 8.5, height = 4.8, device = cairo_pdf)
ggsave(file.path(args$outdir, "panel_B_model_comparison.png"), pB_individual, width = 8.5, height = 4.8, dpi = 450)

# Panel C: GAM improvement heatmap.
gam_sum$taxon_label <- unname(taxon_label_map[gam_sum$taxon])
gam_sum$taxon_label[is.na(gam_sum$taxon_label)] <- gam_sum$taxon[is.na(gam_sum$taxon_label)]
heatmap_order <- taxon_order[taxon_order %in% gam_sum$taxon]
# Reverse factor levels so alphabetical order is read from top to bottom.
gam_sum$taxon_label <- factor(
  gam_sum$taxon_label,
  levels = rev(unname(taxon_label_map[heatmap_order]))
)
pC <- ggplot(gam_sum, aes(x = tool, y = taxon_label, fill = log_error_improvement)) +
  geom_tile() +
  geom_text(
    aes(
      label = calibratability_code,
      colour = ifelse(
        is.na(log_error_improvement),
        "black",
        ifelse(abs(log_error_improvement) >= 0.18, "white", "black")
      )
    ),
    size = 3
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    midpoint = 0,
    name = "Reduction in median absolute\nlog10 error versus raw",
    na.value = "grey90"
  ) +
  labs(
    title = "C. Evidence for profiler–taxon calibratability",
    subtitle = paste0(
      "S = strong; C = conditional; L = little evidence; ",
      "ND = insufficient detection or modelling data"
    ),
    x = NULL,
    y = NULL
  ) + theme_manuscript

pC_individual <- pC + labs(title = NULL, subtitle = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())
ggsave(file.path(args$outdir, "panel_C_calibratability_heatmap.pdf"), pC_individual, width = 7.5, height = 6.5, device = cairo_pdf)
ggsave(file.path(args$outdir, "panel_C_calibratability_heatmap.png"), pC_individual, width = 7.5, height = 6.5, dpi = 450)

# Panel D: cross-study transfer if available.
pD <- NULL
if (nrow(transfer)) {
  trans_agg <- aggregate(transfer[c("median_abs_log_error")],
                         by = transfer[c("tool", "taxon", "model", "train_study", "test_study")], median, na.rm = TRUE)
  trans_agg$model <- factor(trans_agg$model, levels = c("Raw identity", "GAM"))
  clean_study_label <- function(x) gsub("_", " ", x, fixed = TRUE)
  trans_agg$transfer_label <- paste(
    clean_study_label(trans_agg$train_study),
    "→",
    clean_study_label(trans_agg$test_study)
  )
  pD <- ggplot(
    trans_agg,
    aes(
      x = model,
      y = median_abs_log_error,
      group = interaction(tool, taxon, train_study, test_study),
      colour = tool
    )
  ) +
    geom_line(alpha = 0.45, linewidth = 0.45) +
    geom_point(size = 1.7, alpha = 0.9) +
    scale_colour_manual(
      values = c(
        "Kraken2 + Bracken" = "#009E73",
        "MetaPhlAn 4" = "#6F5BD3"
      ),
      name = "Profiler"
    ) +
    facet_wrap(~ transfer_label, scales = "free_y") +
    labs(title = "D. Cross-study transfer",
         subtitle = "Models were trained in one cohort and evaluated in the other",
         x = NULL, y = "Median absolute log10 error") + theme_manuscript
  pD_individual <- pD + labs(title = NULL, subtitle = NULL) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  ggsave(file.path(args$outdir, "panel_D_cross_study_transfer.pdf"), pD_individual, width = 9, height = 4.8, device = cairo_pdf)
  ggsave(file.path(args$outdir, "panel_D_cross_study_transfer.png"), pD_individual, width = 9, height = 4.8, dpi = 450)
}

if (requireNamespace("patchwork", quietly = TRUE)) {
  if (!is.null(pD)) {
    composite <- (pA / (pB | pC) / pD) + patchwork::plot_layout(heights = c(1.25, 1, 0.8))
    h <- 16
  } else {
    composite <- pA / (pB | pC) + patchwork::plot_layout(heights = c(1.3, 1))
    h <- 12
  }
  ggsave(file.path(args$outdir, "manuscript_spike_calibratability_gam.pdf"), composite, width = 14, height = h, device = cairo_pdf)
  ggsave(file.path(args$outdir, "manuscript_spike_calibratability_gam.png"), composite, width = 14, height = h, dpi = 450)
}

writeLines(c(
  paste0("Input: ", args$input),
  paste0("Pseudocount: ", format(pc, scientific = TRUE)),
  paste0("Grouped CV folds: ", args$folds),
  paste0("Grouped CV repeats: ", args$repeats),
  paste0("Minimum overall detection rate for assessment: ", args$min_detection_rate),
  paste0("Minimum detected observations for assessment: ", args$min_detected),
  paste0("Minimum positive-abundance observations for assessment: ", args$min_rows),
  paste0("Minimum biological samples with positive observations: ", effective_min_base_samples),
  paste0("Strong-evidence minimum log-error improvement: > ", args$strong_improvement),
  paste0("Strong-evidence maximum residual median absolute log10 error: < ", args$strong_max_log_error),
  "Primary quantitative models use positive detections only.",
  "Non-detection is summarised separately and is not treated as a calibratable positive abundance.",
  "ND denotes insufficient detection, insufficient modelling data, or failure to obtain a valid grouped cross-validation estimate.",
  "S requires improvement above the strong-improvement threshold and residual error below the strong-error threshold.",
  "C denotes positive held-out improvement that does not satisfy both strong-evidence criteria.",
  "L denotes no positive held-out improvement over the raw identity prediction.",
  "The inverse models predict expected log abundance from observed log abundance.",
  "Figure labels use the short paper label codes (for example Bfrag, Fnuc, Pmic) and are ordered alphabetically by those codes. Panel A displays original-scale abundances on logarithmic axes, with GAMs fitted on the log10 scale and back-transformed for plotting.",
  "Interpretation: proof of calibratability, not a validated correction of natural cohort abundances."
), file.path(args$outdir, "analysis_notes.txt"))

message("[DONE] Outputs written to: ", args$outdir)
