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
cohort_arg <- value("--cohort"); population_arg <- value("--population")
arm_arg <- value("--assembly-arm", "original")
if (any(vapply(list(input, outdir, cohort_arg, population_arg), is.null, logical(1)))) {
  stop("Required: --input FILE --outdir DIR --cohort NAME --population NAME [--assembly-arm NAME]")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dat <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("cohort", "sample_id", "condition", "analysis_population",
              "target_label", "assembly_arm", "profiler",
              "spike_fraction_target", "recovered_spike_signal")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing endpoint columns: ", paste(missing, collapse = ", "))
dat <- dat[dat$cohort == cohort_arg & dat$analysis_population == population_arg &
             dat$assembly_arm == arm_arg, , drop = FALSE]
if (!nrow(dat)) stop("No rows match requested cohort/population/assembly arm.")
dat$dose_value <- as.numeric(dat$spike_fraction_target)
dat$recovered <- as.numeric(dat$recovered_spike_signal)
if (anyNA(dat$dose_value) || any(dat$dose_value <= 0)) stop("Positive finite doses required.")
if (anyNA(dat$recovered) || any(!is.finite(dat$recovered))) stop("Recovered signal must be finite.")
if (length(unique(dat$profiler)) != 2L) stop("Exactly two profilers are required.")
if (length(unique(dat$sample_id)) < 2L) stop("At least two samples are required.")
if (length(unique(dat$target_label)) < 2L) stop("At least two targets are required.")
if (length(unique(dat$dose_value)) < 3L) stop("At least three positive doses are required.")

fmt_dose <- function(x) format(x, scientific = TRUE, digits = 12, trim = TRUE)
dat$dose_pp <- 100 * dat$dose_value
dat$dose_factor <- factor(fmt_dose(dat$dose_value),
                          levels = fmt_dose(sort(unique(dat$dose_value))))
dat$profiler <- factor(dat$profiler, levels = c("kraken2_bracken", "metaphlan4"))
dat$condition <- factor(dat$condition); dat$sample_id <- factor(dat$sample_id)
dat$target_label <- factor(dat$target_label)
linear_formula <- recovered ~ profiler * dose_pp + condition +
  s(sample_id, bs = "re") + s(target_label, bs = "re")
categorical_formula <- recovered ~ profiler * dose_factor + condition +
  s(sample_id, bs = "re") + s(target_label, bs = "re")
linear <- mgcv::gam(linear_formula, data = dat, family = gaussian(), method = "ML")
categorical <- mgcv::gam(categorical_formula, data = dat, family = gaussian(), method = "ML")

model_ok <- function(model) {
  full_convergence <- isTRUE(model$outer.info$conv == "full convergence")
  hessian <- model$outer.info$hess
  hessian_pd <- is.matrix(hessian) && all(is.finite(hessian)) &&
    all(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  finite <- all(is.finite(coef(model)))
  list(converged = isTRUE(model$converged), full_convergence = full_convergence,
       finite_coefficients = finite, hessian_positive_definite = hessian_pd,
       pass = isTRUE(model$converged) && full_convergence && finite && hessian_pd)
}
linear_ok <- model_ok(linear); categorical_ok <- model_ok(categorical)
diagnostics <- rbind(
  data.frame(model = "linear", metric = names(linear_ok)[1:4], value = unlist(linear_ok[1:4])),
  data.frame(model = "categorical", metric = names(categorical_ok)[1:4], value = unlist(categorical_ok[1:4])),
  data.frame(model = "data", metric = c("rows", "samples", "targets", "positive_doses"),
             value = c(nrow(dat), nlevels(dat$sample_id), nlevels(dat$target_label), nlevels(dat$dose_factor)))
)
write.table(diagnostics, file.path(outdir, "continuous_model_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!linear_ok$pass || !categorical_ok$pass) stop("Continuous model failed convergence diagnostics.")

capture.output(summary(linear), file = file.path(outdir, "continuous_linear_model_summary.txt"))
coef_table <- as.data.frame(summary(linear)$p.table)
coef_table$term <- rownames(coef_table); rownames(coef_table) <- NULL
write.table(coef_table, file.path(outdir, "continuous_linear_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
nonlinear_test <- anova(linear, categorical, test = "Chisq")
capture.output(nonlinear_test, file = file.path(outdir, "categorical_dose_nonlinearity_test.txt"))
nonlinear_table <- as.data.frame(nonlinear_test)
nonlinear_table$model <- rownames(nonlinear_table); rownames(nonlinear_table) <- NULL
write.table(nonlinear_table, file.path(outdir, "categorical_dose_nonlinearity_test.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

beta <- coef(linear); covariance <- vcov(linear); base_term <- "dose_pp"
interaction_term <- grep("^profilermetaphlan4:dose_pp$|^dose_pp:profilermetaphlan4$",
                         names(beta), value = TRUE)
if (!base_term %in% names(beta) || length(interaction_term) != 1L) stop("Cannot identify slopes.")
contrast_vectors <- list(
  kraken2_bracken = setNames(as.numeric(names(beta) == base_term), names(beta)),
  metaphlan4 = setNames(as.numeric(names(beta) == base_term |
                                     names(beta) == interaction_term), names(beta)))
slopes <- do.call(rbind, lapply(names(contrast_vectors), function(profiler_i) {
  contrast <- contrast_vectors[[profiler_i]]
  estimate_pp <- sum(contrast * beta)
  se_pp <- sqrt(drop(contrast %*% covariance %*% contrast))
  response_slope <- estimate_pp / 0.01; response_se <- se_pp / 0.01
  z <- (response_slope - 1) / response_se
  data.frame(profiler = profiler_i, slope_per_percentage_point = estimate_pp,
             slope_se_per_percentage_point = se_pp, response_slope = response_slope,
             response_slope_se = response_se,
             lower_95 = response_slope - 1.96 * response_se,
             upper_95 = response_slope + 1.96 * response_se,
             deviation_from_one = response_slope - 1,
             p_value_vs_one = 2 * pnorm(abs(z), lower.tail = FALSE))
}))
slopes$q_value_bh_vs_one <- p.adjust(slopes$p_value_vs_one, method = "BH")
write.table(slopes, file.path(outdir, "profiler_response_slopes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

interaction_name <- interaction_term[[1]]
interaction <- data.frame(
  term = interaction_name,
  difference_response_slope_metaphlan_minus_bracken = beta[[interaction_name]] / 0.01,
  standard_error = sqrt(covariance[interaction_name, interaction_name]) / 0.01,
  p_value = coef_table[match(interaction_name, coef_table$term), "Pr(>|t|)"])
interaction$lower_95 <- interaction$difference_response_slope_metaphlan_minus_bracken - 1.96 * interaction$standard_error
interaction$upper_95 <- interaction$difference_response_slope_metaphlan_minus_bracken + 1.96 * interaction$standard_error
write.table(interaction, file.path(outdir, "primary_profiler_slope_contrast.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

residual_rows <- data.frame(
  cohort = cohort_arg, population = population_arg, assembly_arm = arm_arg,
  sample_id = dat$sample_id, target_label = dat$target_label, profiler = dat$profiler,
  spike_fraction_target = dat$dose_value, recovered_spike_signal = dat$recovered,
  fitted = fitted(linear), residual = residuals(linear, type = "response"))
write.table(residual_rows, file.path(outdir, "continuous_fitted_residuals.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
groups <- split(residual_rows, list(residual_rows$profiler,
                                    residual_rows$spike_fraction_target), drop = TRUE)
residual_summary <- do.call(rbind, lapply(groups, function(z) data.frame(
  profiler = as.character(z$profiler[1]), spike_fraction_target = z$spike_fraction_target[1],
  n = nrow(z), residual_mean = mean(z$residual), residual_sd = sd(z$residual),
  residual_median = median(z$residual), residual_mad = mad(z$residual, constant = 1))))
write.table(residual_summary, file.path(outdir, "continuous_residuals_by_dose.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
seal_files <- c(normalizePath(input), file.path(outdir, "continuous_model_diagnostics.tsv"),
  file.path(outdir, "continuous_linear_coefficients.tsv"),
  file.path(outdir, "profiler_response_slopes.tsv"),
  file.path(outdir, "primary_profiler_slope_contrast.tsv"),
  file.path(outdir, "categorical_dose_nonlinearity_test.tsv"),
  file.path(outdir, "continuous_fitted_residuals.tsv"),
  file.path(outdir, "continuous_residuals_by_dose.tsv"), file.path(outdir, "session_info.txt"))
status <- system2("sha256sum", seal_files, stdout = file.path(outdir, "continuous_model.sha256"))
if (!identical(status, 0L)) stop("Could not create SHA-256 seal.")
writeLines(c(paste("cohort", cohort_arg, sep = "\t"),
  paste("population", population_arg, sep = "\t"),
  paste("assembly_arm", arm_arg, sep = "\t"), paste("rows", nrow(dat), sep = "\t"),
  "status\tPASS"), file.path(outdir, "SUCCESS"))
message("[PASS] Continuous dose-response model completed: ", outdir)
