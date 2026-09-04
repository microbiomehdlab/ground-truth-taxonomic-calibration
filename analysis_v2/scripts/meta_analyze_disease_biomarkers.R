#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args); if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
input_paths <- strsplit(value("--inputs", ""), ",", fixed = TRUE)[[1]]
expected <- strsplit(value("--expected-cohorts", "yachida,feng,zeller"), ",", fixed = TRUE)[[1]]
outdir <- value("--outdir")
if (!length(input_paths) || any(!nzchar(input_paths)) || is.null(outdir))
  stop("Required: --inputs FILE[,FILE...] --outdir DIR")
if (length(expected) < 2L || anyDuplicated(expected)) stop("Expected cohorts must be distinct.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

needed <- c("cohort", "study", "analysis_population", "target_label", "assembly_arm",
            "profiler", "dose_level", "spike_fraction_target", "contrast", "feature",
            "effect", "standard_error", "p_value", "q_value", "model_spec", "include")
tables <- lapply(input_paths, function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(needed, names(x)); if (length(missing)) stop("Missing columns in ", path)
  x$source_result <- normalizePath(path); x
})
data <- do.call(rbind, tables)
data <- data[data$include == 1 & data$model_spec == "primary_age_sex", , drop = FALSE]
data$effect <- as.numeric(data$effect); data$standard_error <- as.numeric(data$standard_error)
data$spike_fraction_target <- as.numeric(data$spike_fraction_target)
if (!nrow(data) || anyNA(data$effect) || anyNA(data$standard_error) ||
    any(data$standard_error < 0)) stop("Invalid cohort estimates.")
observed <- sort(unique(data$cohort))
if (!setequal(observed, expected))
  stop("Observed cohorts do not equal expected cohorts: observed=", paste(observed, collapse = ","),
       " expected=", paste(sort(expected), collapse = ","))

key <- c("analysis_population", "target_label", "assembly_arm", "profiler",
         "dose_level", "contrast", "feature", "model_spec")
identity <- do.call(paste, c(data[c("cohort", key)], sep = "\r"))
if (anyDuplicated(identity)) stop("Duplicate cohort estimate for a synthesis cell.")

reml <- function(y, se) {
  v <- se^2; k <- length(y)
  objective <- function(tau2) {
    w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
    sum(log(v + tau2)) + log(sum(w)) + sum(w * (y - mu)^2)
  }
  upper <- max(var(y, na.rm = TRUE) * 100, max(v) * 100, 1)
  tau2 <- optimize(objective, c(0, upper))$minimum
  if (tau2 < upper * 1e-8) tau2 <- 0
  w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
  q_re <- sum(w * (y - mu)^2) / (k - 1)
  se_hk <- sqrt(max(1, q_re) / sum(w))
  crit <- qt(.975, k - 1); statistic <- mu / se_hk
  p <- max(2 * pt(abs(statistic), k - 1, lower.tail = FALSE), .Machine$double.xmin)
  q <- sum((1 / v) * (y - sum(y / v) / sum(1 / v))^2)
  i2 <- if (q > 0) max(0, (q - (k - 1)) / q) * 100 else 0
  prediction_crit <- if (k >= 3L) qt(.975, k - 2) else NA_real_
  prediction_half <- prediction_crit * sqrt(tau2 + se_hk^2)
  c(effect = mu, standard_error = se_hk, lower_95 = mu - crit * se_hk,
    upper_95 = mu + crit * se_hk, p_value = p, tau2 = tau2, Q = q,
    Q_df = k - 1, I2_pct = i2,
    prediction_lower_95 = mu - prediction_half,
    prediction_upper_95 = mu + prediction_half)
}

groups <- split(data, interaction(data[key], drop = TRUE, lex.order = TRUE))
pooled <- list(); coverage <- list(); loo <- list(); pi <- ci <- li <- 1L
for (group in groups) {
  first <- group[1, , drop = FALSE]
  present <- sort(unique(group$cohort)); missing <- setdiff(expected, present)
  usable <- is.finite(group$standard_error) & group$standard_error > 0
  status <- if (length(missing)) "INCOMPLETE_COHORT_COVERAGE" else if (!all(usable)) "NONPOSITIVE_STANDARD_ERROR" else "ELIGIBLE"
  coverage[[ci]] <- c(as.list(first[1, key, drop = FALSE]), list(
    cohorts_present = paste(present, collapse = ","), cohorts_missing = paste(missing, collapse = ","),
    cohort_count = length(present), status = status)); ci <- ci + 1L
  if (status != "ELIGIBLE") next
  stats <- reml(group$effect, group$standard_error)
  signs <- sign(group$effect); nonzero <- signs[signs != 0]
  pooled[[pi]] <- c(as.list(first[1, key, drop = FALSE]), as.list(stats), list(
    cohort_count = nrow(group), cohorts = paste(sort(group$cohort), collapse = ","),
    dose_fraction_median = median(group$spike_fraction_target),
    dose_fraction_min = min(group$spike_fraction_target), dose_fraction_max = max(group$spike_fraction_target),
    direction_consistent = as.integer(length(unique(nonzero)) <= 1L), q_value = NA_real_)); pi <- pi + 1L
  for (held in sort(group$cohort)) {
    retained <- group[group$cohort != held, , drop = FALSE]
    fixed_w <- 1 / retained$standard_error^2
    estimate <- sum(fixed_w * retained$effect) / sum(fixed_w)
    loo[[li]] <- c(as.list(first[1, key, drop = FALSE]), list(
      omitted_cohort = held, retained_cohorts = paste(sort(retained$cohort), collapse = ","),
      fixed_effect_estimate = estimate)); li <- li + 1L
  }
}
bind <- function(x) if (length(x)) do.call(rbind.data.frame, c(x, stringsAsFactors = FALSE)) else data.frame()
pooled <- bind(pooled); coverage <- bind(coverage); loo <- bind(loo)
if (!nrow(pooled)) stop("No complete cross-cohort synthesis cells.")
bh_key <- c("analysis_population", "target_label", "assembly_arm", "profiler",
            "dose_level", "contrast", "model_spec")
families <- split(seq_len(nrow(pooled)), interaction(pooled[bh_key], drop = TRUE, lex.order = TRUE))
for (indices in families) pooled$q_value[indices] <- p.adjust(as.numeric(pooled$p_value[indices]), "BH")
write.table(data, file.path(outdir, "cohort_specific_disease_estimates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pooled, file.path(outdir, "random_effects_meta_analysis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coverage, file.path(outdir, "synthesis_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(loo, file.path(outdir, "leave_one_cohort_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
settings <- data.frame(setting = c("estimator", "uncertainty", "prediction_interval", "minimum_cohorts",
                                   "eligibility", "multiplicity"),
  value = c("REML random effects", "modified Hartung-Knapp (scale >= 1)",
            "t interval with k-2 df", length(expected), "present with positive SE in every expected cohort",
            "BH across features within population/target/arm/profiler/dose/contrast/model"))
write.table(settings, file.path(outdir, "meta_analysis_settings.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
summary <- data.frame(metric = c("expected_cohorts", "eligible_cells", "coverage_rows", "status"),
                      value = c(length(expected), nrow(pooled), nrow(coverage), "PASS"))
write.table(summary, file.path(outdir, "meta_analysis_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
files <- c(normalizePath(input_paths), list.files(outdir, full.names = TRUE))
status <- system2("sha256sum", files, stdout = file.path(outdir, "meta_analysis.sha256"))
if (!identical(status, 0L)) stop("Could not seal meta-analysis.")
writeLines(c("analysis\tcross_cohort_disease_biomarker_meta_analysis", "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Cross-cohort meta-analysis completed: ", outdir)
