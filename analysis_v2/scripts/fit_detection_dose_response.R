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

input <- value("--input")
outdir <- value("--outdir")
cohort_arg <- value("--cohort")
population_arg <- value("--population")
arm_arg <- value("--assembly-arm", "original")
if (is.null(input) || is.null(outdir) || is.null(cohort_arg) || is.null(population_arg)) {
  stop("Required: --input FILE --outdir DIR --cohort NAME --population NAME [--assembly-arm NAME]")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dat <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
required <- c(
  "cohort", "sample_id", "condition", "analysis_population", "target_label",
  "assembly_arm", "profiler", "spike_fraction_target"
)
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing endpoint columns: ", paste(missing, collapse = ", "))

dat <- dat[
  dat$cohort == cohort_arg &
  dat$analysis_population == population_arg &
  dat$assembly_arm == arm_arg,
  , drop = FALSE
]
if ("include" %in% names(dat)) dat <- dat[dat$include == 1, , drop = FALSE]
if (!nrow(dat)) stop("No rows match requested cohort/population/assembly arm.")
detection_field <- if ("detected_native_nonzero" %in% names(dat)) {
  "detected_native_nonzero"
} else if ("observed_detected_native_nonzero" %in% names(dat)) {
  "observed_detected_native_nonzero"
} else {
  stop("Missing detected_native_nonzero column.")
}
dat$detected <- as.integer(dat[[detection_field]])
dat$dose_value <- as.numeric(dat$spike_fraction_target)
if (anyNA(dat$detected) || any(!dat$detected %in% 0:1)) stop("Detection must be binary.")
if (anyNA(dat$dose_value) || any(dat$dose_value < 0)) stop("Invalid dose values.")

dose_values <- sort(unique(dat$dose_value))
if (!0 %in% dose_values || length(dose_values) < 2L) stop("At least zero and one positive dose are required.")
if (length(unique(dat$profiler)) != 2L) stop("Exactly two profilers are required.")
if (length(unique(dat$sample_id)) < 2L) stop("At least two biological samples are required.")
if (length(unique(dat$target_label)) < 2L) stop("At least two targets are required.")
if (length(unique(dat$detected)) != 2L) stop("Both detection outcome classes are required.")

fmt_dose <- function(x) format(x, scientific = TRUE, digits = 12, trim = TRUE)
dat$dose <- factor(fmt_dose(dat$dose_value), levels = fmt_dose(dose_values))
dat$profiler <- factor(dat$profiler, levels = c("kraken2_bracken", "metaphlan4"))
dat$condition <- factor(dat$condition)
dat$sample_id <- factor(dat$sample_id)
dat$target_label <- factor(dat$target_label)

full_formula <- detected ~ profiler * dose + condition +
  s(sample_id, bs = "re") + s(target_label, bs = "re")
reduced_formula <- detected ~ profiler + dose + condition +
  s(sample_id, bs = "re") + s(target_label, bs = "re")
full <- mgcv::gam(full_formula, data = dat, family = binomial(), method = "ML")
reduced <- mgcv::gam(reduced_formula, data = dat, family = binomial(), method = "ML")

full_convergence <- isTRUE(full$outer.info$conv == "full convergence")
hessian <- full$outer.info$hess
hessian_pd <- is.matrix(hessian) && all(is.finite(hessian)) &&
  all(eigen(hessian, symmetric = TRUE, only.values = TRUE)$values > 1e-8)

diagnostic <- data.frame(
  metric = c("rows", "samples", "targets", "doses", "converged", "full_convergence",
             "finite_coefficients", "hessian_positive_definite"),
  value = c(nrow(dat), length(levels(dat$sample_id)), length(levels(dat$target_label)),
            length(levels(dat$dose)), isTRUE(full$converged), full_convergence,
            all(is.finite(coef(full))), hessian_pd),
  stringsAsFactors = FALSE
)
write.table(diagnostic, file.path(outdir, "detection_model_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!isTRUE(full$converged) || !full_convergence || any(!is.finite(coef(full))) ||
    !hessian_pd) {
  stop("Primary detection model failed convergence diagnostics.")
}

comparison <- anova(reduced, full, test = "Chisq")
capture.output(summary(full), file = file.path(outdir, "detection_model_summary.txt"))
capture.output(comparison, file = file.path(outdir, "profiler_by_dose_primary_test.txt"))
comparison_table <- as.data.frame(comparison)
comparison_table$model <- rownames(comparison_table)
rownames(comparison_table) <- NULL
write.table(comparison_table, file.path(outdir, "profiler_by_dose_primary_test.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
coef_table <- as.data.frame(summary(full)$p.table)
coef_table$term <- rownames(coef_table)
rownames(coef_table) <- NULL
write.table(coef_table, file.path(outdir, "detection_model_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Average fixed-effect predictions over the observed condition distribution.
conditions <- levels(dat$condition)
weights <- prop.table(table(dat$condition))[conditions]
grid <- expand.grid(
  profiler = levels(dat$profiler), dose = levels(dat$dose), condition = conditions,
  stringsAsFactors = FALSE
)
grid$profiler <- factor(grid$profiler, levels = levels(dat$profiler))
grid$dose <- factor(grid$dose, levels = levels(dat$dose))
grid$condition <- factor(grid$condition, levels = levels(dat$condition))
grid$sample_id <- factor(levels(dat$sample_id)[1], levels = levels(dat$sample_id))
grid$target_label <- factor(levels(dat$target_label)[1], levels = levels(dat$target_label))
random_terms <- c("s(sample_id)", "s(target_label)")
link <- predict(full, newdata = grid, type = "link", se.fit = TRUE, exclude = random_terms)
grid$fit_link <- as.numeric(link$fit)
grid$se_link <- as.numeric(link$se.fit)
grid$probability <- plogis(grid$fit_link)
grid$lower_95 <- plogis(grid$fit_link - 1.96 * grid$se_link)
grid$upper_95 <- plogis(grid$fit_link + 1.96 * grid$se_link)
grid$weight <- as.numeric(weights[as.character(grid$condition)])

weighted <- function(x, w) sum(x * w) / sum(w)
curve_rows <- do.call(rbind, lapply(split(grid, list(grid$profiler, grid$dose), drop = TRUE), function(z) {
  data.frame(
    profiler = as.character(z$profiler[1]),
    dose = as.character(z$dose[1]),
    spike_fraction_target = dose_values[match(as.character(z$dose[1]), fmt_dose(dose_values))],
    predicted_detection_probability = weighted(z$probability, z$weight),
    lower_95_descriptive = weighted(z$lower_95, z$weight),
    upper_95_descriptive = weighted(z$upper_95, z$weight)
  )
}))
write.table(curve_rows, file.path(outdir, "detection_probability_curve.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Fixed-effect link-scale profiler contrasts at each dose, condition-weighted.
xp <- predict(full, newdata = grid, type = "lpmatrix", exclude = random_terms)
contrasts <- lapply(levels(dat$dose), function(dose_i) {
  a <- grid$dose == dose_i & grid$profiler == "kraken2_bracken"
  b <- grid$dose == dose_i & grid$profiler == "metaphlan4"
  xa <- colSums(xp[a, , drop = FALSE] * grid$weight[a]) / sum(grid$weight[a])
  xb <- colSums(xp[b, , drop = FALSE] * grid$weight[b]) / sum(grid$weight[b])
  delta <- xb - xa
  estimate <- sum(delta * coef(full))
  se <- sqrt(drop(delta %*% vcov(full) %*% delta))
  z <- estimate / se
  data.frame(
    dose = dose_i,
    spike_fraction_target = dose_values[match(dose_i, fmt_dose(dose_values))],
    log_odds_ratio_metaphlan_vs_bracken = estimate,
    standard_error = se,
    odds_ratio = exp(estimate),
    lower_95 = exp(estimate - 1.96 * se),
    upper_95 = exp(estimate + 1.96 * se),
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
})
contrasts <- do.call(rbind, contrasts)
positive <- contrasts$spike_fraction_target > 0
contrasts$q_value_bh_positive_doses <- NA_real_
contrasts$q_value_bh_positive_doses[positive] <- p.adjust(contrasts$p_value[positive], method = "BH")
write.table(contrasts, file.path(outdir, "profiler_contrasts_by_dose.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
seal_files <- c(
  normalizePath(input),
  file.path(outdir, "detection_model_diagnostics.tsv"),
  file.path(outdir, "detection_model_coefficients.tsv"),
  file.path(outdir, "detection_probability_curve.tsv"),
  file.path(outdir, "profiler_contrasts_by_dose.tsv"),
  file.path(outdir, "profiler_by_dose_primary_test.tsv"),
  file.path(outdir, "session_info.txt")
)
checksum_status <- system2("sha256sum", seal_files,
                           stdout = file.path(outdir, "detection_model.sha256"))
if (!identical(checksum_status, 0L)) stop("Could not create SHA-256 seal.")

writeLines(c(
  paste("cohort", cohort_arg, sep = "\t"),
  paste("population", population_arg, sep = "\t"),
  paste("assembly_arm", arm_arg, sep = "\t"),
  paste("rows", nrow(dat), sep = "\t"),
  "status\tPASS"
), file.path(outdir, "SUCCESS"))
message("[PASS] Detection dose-response model completed: ", outdir)
