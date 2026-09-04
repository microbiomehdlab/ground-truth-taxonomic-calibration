#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' is required.")
})
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
input <- value("--input"); outdir <- value("--outdir")
cohort_arg <- value("--cohort", "yachida")
population_arg <- value("--population", "independent")
if (is.null(input) || is.null(outdir)) {
  stop("Required: --input FILE --outdir DIR [--cohort NAME] [--population NAME]")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dat <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("cohort", "sample_id", "condition", "analysis_population",
              "target_label", "assembly_arm", "profiler",
              "spike_fraction_target", "recovered_spike_signal")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing endpoint columns: ", paste(missing, collapse = ", "))
dat <- dat[dat$cohort == cohort_arg & dat$analysis_population == population_arg &
             dat$assembly_arm %in% c("original", "clean"), , drop = FALSE]
if (!nrow(dat)) stop("No rows match the requested assembly-sensitivity analysis.")
dat$dose_value <- as.numeric(dat$spike_fraction_target)
dat$recovered <- as.numeric(dat$recovered_spike_signal)
if (anyNA(dat$dose_value) || any(dat$dose_value <= 0)) stop("Positive finite doses required.")
if (anyNA(dat$recovered) || any(!is.finite(dat$recovered))) stop("Recovered signal must be finite.")
if (!identical(sort(unique(dat$profiler)), sort(c("kraken2_bracken", "metaphlan4"))))
  stop("Exactly the two frozen profilers are required.")
if (!identical(sort(unique(dat$assembly_arm)), c("clean", "original")))
  stop("Both original and clean assembly arms are required.")
if (length(unique(dat$sample_id)) < 2L || length(unique(dat$target_label)) != 2L)
  stop("At least two samples and exactly two assembly-sensitivity targets are required.")
cells <- xtabs(~ sample_id + target_label + profiler + assembly_arm, data = dat)
if (any(cells == 0L) || length(unique(as.vector(cells))) != 1L)
  stop("The profiler-by-arm design is incomplete or unbalanced.")

dat$dose_pp <- 100 * dat$dose_value
dat$profiler <- factor(dat$profiler, levels = c("kraken2_bracken", "metaphlan4"))
dat$assembly_arm <- factor(dat$assembly_arm, levels = c("original", "clean"))
dat$condition <- factor(dat$condition); dat$sample_id <- factor(dat$sample_id)
dat$target_label <- factor(dat$target_label)
dat$sample_target_profiler <- interaction(
  dat$sample_id, dat$target_label, dat$profiler, drop = TRUE)
dat$clean_dose_pp <- dat$dose_pp * as.numeric(dat$assembly_arm == "clean")
formula <- recovered ~ profiler * assembly_arm * dose_pp * target_label + condition +
  s(sample_id, bs = "re") +
  s(sample_target_profiler, bs = "re") +
  s(sample_target_profiler, by = dose_pp, bs = "re") +
  s(sample_target_profiler, by = clean_dose_pp, bs = "re")
model <- mgcv::gam(formula, data = dat, family = gaussian(), method = "ML")
hessian <- model$outer.info$hess
diagnostics <- data.frame(
  metric = c("converged", "full_convergence", "finite_coefficients",
             "hessian_positive_definite", "rows", "samples", "targets",
             "sample_target_profiler_clusters", "achieved_doses"),
  value = c(isTRUE(model$converged), isTRUE(model$outer.info$conv == "full convergence"),
            all(is.finite(coef(model))), is.matrix(hessian) && all(is.finite(hessian)) &&
              all(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values > 1e-8),
            nrow(dat), nlevels(dat$sample_id), nlevels(dat$target_label),
            nlevels(dat$sample_target_profiler),
            length(unique(dat$dose_value))))
write.table(diagnostics, file.path(outdir, "assembly_model_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(as.logical(diagnostics$value[1:4]))) stop("Assembly model failed convergence diagnostics.")

beta <- coef(model); covariance <- vcov(model)
reference <- data.frame(condition = dat$condition[1], sample_id = dat$sample_id[1],
                        stringsAsFactors = FALSE)
smooth_terms <- c("s(sample_id)", "s(sample_target_profiler)",
                  "s(sample_target_profiler):dose_pp",
                  "s(sample_target_profiler):clean_dose_pp")
slope_contrast <- function(profiler_i, arm_i, target_i) {
  low <- transform(reference, profiler = profiler_i, assembly_arm = arm_i,
                   target_label = target_i, dose_pp = 0, clean_dose_pp = 0)
  high <- transform(reference, profiler = profiler_i, assembly_arm = arm_i,
                    target_label = target_i, dose_pp = 1,
                    clean_dose_pp = as.numeric(arm_i == "clean"))
  low$profiler <- factor(low$profiler, levels = levels(dat$profiler))
  high$profiler <- factor(high$profiler, levels = levels(dat$profiler))
  low$assembly_arm <- factor(low$assembly_arm, levels = levels(dat$assembly_arm))
  high$assembly_arm <- factor(high$assembly_arm, levels = levels(dat$assembly_arm))
  low$condition <- factor(low$condition, levels = levels(dat$condition))
  high$condition <- factor(high$condition, levels = levels(dat$condition))
  low$sample_id <- factor(low$sample_id, levels = levels(dat$sample_id))
  high$sample_id <- factor(high$sample_id, levels = levels(dat$sample_id))
  low$target_label <- factor(low$target_label, levels = levels(dat$target_label))
  high$target_label <- factor(high$target_label, levels = levels(dat$target_label))
  low$sample_target_profiler <- interaction(
    low$sample_id, low$target_label, low$profiler, drop = TRUE)
  high$sample_target_profiler <- interaction(
    high$sample_id, high$target_label, high$profiler, drop = TRUE)
  low$sample_target_profiler <- factor(
    low$sample_target_profiler, levels = levels(dat$sample_target_profiler))
  high$sample_target_profiler <- factor(
    high$sample_target_profiler, levels = levels(dat$sample_target_profiler))
  drop(predict(model, high, type = "lpmatrix", exclude = smooth_terms) -
         predict(model, low, type = "lpmatrix", exclude = smooth_terms)) / 0.01
}
vectors <- list()
for (profiler_i in levels(dat$profiler)) for (arm_i in levels(dat$assembly_arm))
  for (target_i in levels(dat$target_label))
    vectors[[paste(profiler_i, arm_i, target_i, sep = "__")]] <-
      slope_contrast(profiler_i, arm_i, target_i)
estimate <- function(vector) {
  est <- sum(vector * beta); se <- sqrt(drop(vector %*% covariance %*% vector))
  c(estimate = est, standard_error = se, lower_95 = est - 1.96 * se,
    upper_95 = est + 1.96 * se,
    p_value = pmax(2 * pnorm(abs(est / se), lower.tail = FALSE),
                   .Machine$double.xmin))
}
slopes <- do.call(rbind, lapply(names(vectors), function(name) {
  parts <- strsplit(name, "__", fixed = TRUE)[[1]]
  stats <- estimate(vectors[[name]])
  data.frame(profiler = parts[1], assembly_arm = parts[2], target_label = parts[3],
             response_slope = unname(stats[1]),
             standard_error = unname(stats[2]), lower_95 = unname(stats[3]),
             upper_95 = unname(stats[4]), deviation_from_one = unname(stats[1]) - 1,
             p_value_vs_one = pmax(
               2 * pnorm(abs((unname(stats[1]) - 1) / unname(stats[2])),
                         lower.tail = FALSE), .Machine$double.xmin))
}))
slopes$q_value_bh_vs_one <- p.adjust(slopes$p_value_vs_one, method = "BH")
slopes$minus_log10_p_value_vs_one <- -log10(slopes$p_value_vs_one)
write.table(slopes, file.path(outdir, "assembly_profiler_response_slopes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

contrasts <- list()
for (target_i in levels(dat$target_label)) {
  for (profiler_i in levels(dat$profiler)) {
    name <- paste(profiler_i, target_i, "clean_minus_original", sep = "__")
    contrasts[[name]] <- vectors[[paste(profiler_i, "clean", target_i, sep = "__")]] -
      vectors[[paste(profiler_i, "original", target_i, sep = "__")]]
  }
  contrasts[[paste("profiler_difference_in_differences", target_i, sep = "__")]] <-
    contrasts[[paste("metaphlan4", target_i, "clean_minus_original", sep = "__")]] -
    contrasts[[paste("kraken2_bracken", target_i, "clean_minus_original", sep = "__")]]
}
for (profiler_i in levels(dat$profiler)) {
  contrasts[[paste(profiler_i, "pooled_clean_minus_original", sep = "__")]] <-
    Reduce(`+`, lapply(levels(dat$target_label), function(target_i)
      contrasts[[paste(profiler_i, target_i, "clean_minus_original", sep = "__")]])) /
    nlevels(dat$target_label)
}
contrasts[["pooled_profiler_difference_in_differences"]] <-
  contrasts[["metaphlan4__pooled_clean_minus_original"]] -
  contrasts[["kraken2_bracken__pooled_clean_minus_original"]]
contrast_table <- do.call(rbind, lapply(names(contrasts), function(name) {
  stats <- estimate(contrasts[[name]])
  data.frame(contrast = name, estimate = unname(stats[1]), standard_error = unname(stats[2]),
             lower_95 = unname(stats[3]), upper_95 = unname(stats[4]),
             p_value = unname(stats[5]))
}))
contrast_table$contrast_family <- ifelse(
  grepl("^(kraken2_bracken|metaphlan4)__(Pana|Pint)__clean_minus_original$",
        contrast_table$contrast), "primary_target_specific_assembly",
  ifelse(grepl("^profiler_difference_in_differences__(Pana|Pint)$",
               contrast_table$contrast), "secondary_target_specific_profiler_did",
  ifelse(grepl("^(kraken2_bracken|metaphlan4)__pooled_clean_minus_original$",
               contrast_table$contrast), "secondary_pooled_assembly",
         "secondary_pooled_profiler_did")))
contrast_table$q_value_bh <- ave(
  contrast_table$p_value, contrast_table$contrast_family,
  FUN = function(values) p.adjust(values, method = "BH"))
contrast_table$minus_log10_p_value <- -log10(contrast_table$p_value)
write.table(contrast_table, file.path(outdir, "assembly_slope_contrasts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

capture.output(summary(model), file = file.path(outdir, "assembly_model_summary.txt"))
fitted_rows <- data.frame(
  cohort = cohort_arg, sample_id = dat$sample_id, condition = dat$condition,
  target_label = dat$target_label, profiler = dat$profiler, assembly_arm = dat$assembly_arm,
  spike_fraction_target = dat$dose_value, recovered_spike_signal = dat$recovered,
  fitted = fitted(model), residual = residuals(model, type = "response"))
write.table(fitted_rows, file.path(outdir, "assembly_fitted_residuals.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
seal_files <- c(normalizePath(input), file.path(outdir, "assembly_model_diagnostics.tsv"),
                file.path(outdir, "assembly_profiler_response_slopes.tsv"),
                file.path(outdir, "assembly_slope_contrasts.tsv"),
                file.path(outdir, "assembly_fitted_residuals.tsv"), file.path(outdir, "session_info.txt"))
status <- system2("sha256sum", seal_files, stdout = file.path(outdir, "assembly_model.sha256"))
if (!identical(status, 0L)) stop("Could not create SHA-256 seal.")
writeLines(c(paste("cohort", cohort_arg, sep = "\t"), paste("population", population_arg, sep = "\t"),
             paste("rows", nrow(dat), sep = "\t"), "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Assembly-choice sensitivity model completed: ", outdir)
