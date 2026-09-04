#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value for ", flag)
  args[[hit + 1L]]
}
input <- value("--input")
outdir <- value("--outdir")
cohort_arg <- value("--cohort", "yachida")
population_arg <- value("--population", "independent")
bootstrap_replicates <- as.integer(value("--bootstrap-replicates", "20000"))
signflip_replicates <- as.integer(value("--signflip-replicates", "100000"))
seed <- as.integer(value("--seed", "20260904"))
if (is.null(input) || is.null(outdir)) {
  stop("Required: --input FILE --outdir DIR [--cohort NAME] [--population NAME]")
}
if (bootstrap_replicates < 999L || signflip_replicates < 999L || is.na(seed)) {
  stop("At least 999 bootstrap and sign-flip replicates and a finite seed are required.")
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
dat$dose <- as.numeric(dat$spike_fraction_target)
dat$response <- as.numeric(dat$recovered_spike_signal)
if (!nrow(dat) || anyNA(dat$dose) || any(dat$dose <= 0) ||
    anyNA(dat$response) || any(!is.finite(dat$response))) {
  stop("Finite response rows at positive doses are required.")
}
if (!identical(sort(unique(dat$profiler)), sort(c("kraken2_bracken", "metaphlan4"))) ||
    !identical(sort(unique(dat$assembly_arm)), c("clean", "original")) ||
    length(unique(dat$target_label)) != 2L) {
  stop("The frozen two-profiler, two-arm, two-target design is required.")
}

key_columns <- c("sample_id", "condition", "target_label", "profiler", "assembly_arm")
keys <- interaction(dat[key_columns], drop = TRUE, lex.order = TRUE)
groups <- split(dat, keys)
expected_doses <- sort(unique(dat$dose))
if (length(expected_doses) < 3L) stop("At least three positive doses are required.")

fit_one <- function(frame) {
  if (nrow(frame) != length(expected_doses) ||
      !identical(sort(frame$dose), expected_doses) || anyDuplicated(frame$dose)) {
    stop("Every sample-target-profiler-arm cell must contain each dose exactly once.")
  }
  model <- lm(response ~ dose, data = frame)
  first <- frame[1L, key_columns, drop = FALSE]
  data.frame(first, n_doses = nrow(frame), intercept = unname(coef(model)[1L]),
             response_slope = unname(coef(model)[2L]), r_squared = summary(model)$r.squared,
             stringsAsFactors = FALSE)
}
sample_slopes <- do.call(rbind, lapply(groups, fit_one))
row.names(sample_slopes) <- NULL

wide <- reshape(sample_slopes[c("sample_id", "condition", "target_label", "profiler",
                                "assembly_arm", "response_slope")],
                idvar = c("sample_id", "condition", "target_label", "profiler"),
                timevar = "assembly_arm", direction = "wide")
needed <- c("response_slope.original", "response_slope.clean")
if (!all(needed %in% names(wide)) || anyNA(wide[needed])) stop("Arm pairing failed.")
names(wide)[match(needed, names(wide))] <- c("original_slope", "clean_slope")
wide$clean_minus_original <- wide$clean_slope - wide$original_slope
paired_differences <- wide

set.seed(seed)
infer <- function(values) {
  values <- as.numeric(values); n <- length(values)
  if (n < 2L || any(!is.finite(values))) stop("At least two finite sample-level values required.")
  boot <- replicate(bootstrap_replicates, mean(sample(values, n, replace = TRUE)))
  ci <- unname(quantile(boot, c(0.025, 0.975), type = 6))
  observed <- abs(mean(values))
  if (n <= 16L) {
    codes <- 0:(2^n - 1L)
    null <- vapply(codes, function(code) {
      signs <- ifelse(bitwAnd(code, bitwShiftL(1L, 0:(n - 1L))) == 0L, -1, 1)
      mean(signs * values)
    }, numeric(1))
    p <- mean(abs(null) >= observed - 1e-15)
    method <- "exact_sign_flip"
    permutations <- length(null)
  } else {
    exceed <- 0L
    for (i in seq_len(signflip_replicates)) {
      exceed <- exceed + as.integer(abs(mean(sample(c(-1, 1), n, replace = TRUE) * values)) >= observed)
    }
    p <- (exceed + 1) / (signflip_replicates + 1)
    method <- "monte_carlo_sign_flip"
    permutations <- signflip_replicates
  }
  c(n_samples = n, mean = mean(values), median = median(values), sd = sd(values),
    lower_95_bootstrap = ci[1L], upper_95_bootstrap = ci[2L], p_value = p,
    sign_flip_replicates = permutations, inference_method = method)
}

summarize_groups <- function(frame, group_columns, value_column, family) {
  grouping <- split(frame, interaction(frame[group_columns], drop = TRUE, lex.order = TRUE))
  result <- do.call(rbind, lapply(grouping, function(part) {
    stats <- infer(part[[value_column]])
    data.frame(part[1L, group_columns, drop = FALSE],
               contrast_family = family, estimand = value_column,
               n_samples = as.integer(stats[["n_samples"]]),
               mean = as.numeric(stats[["mean"]]), median = as.numeric(stats[["median"]]),
               sd = as.numeric(stats[["sd"]]),
               lower_95_bootstrap = as.numeric(stats[["lower_95_bootstrap"]]),
               upper_95_bootstrap = as.numeric(stats[["upper_95_bootstrap"]]),
               p_value = as.numeric(stats[["p_value"]]),
               inference_method = stats[["inference_method"]],
               sign_flip_replicates = as.integer(stats[["sign_flip_replicates"]]),
               stringsAsFactors = FALSE)
  }))
  row.names(result) <- NULL
  result$q_value_bh <- p.adjust(result$p_value, method = "BH")
  result
}

primary <- summarize_groups(paired_differences, c("target_label", "profiler"),
                            "clean_minus_original", "primary_target_profiler")

did <- reshape(paired_differences[c("sample_id", "condition", "target_label", "profiler",
                                    "clean_minus_original")],
               idvar = c("sample_id", "condition", "target_label"),
               timevar = "profiler", direction = "wide")
did$metaphlan_minus_kraken <- did$clean_minus_original.metaphlan4 -
  did$clean_minus_original.kraken2_bracken
profiler_did <- summarize_groups(did, "target_label", "metaphlan_minus_kraken",
                                 "secondary_target_profiler_did")

pooled <- aggregate(clean_minus_original ~ sample_id + condition + profiler,
                    paired_differences, mean)
pooled_summary <- summarize_groups(pooled, "profiler", "clean_minus_original",
                                   "secondary_equal_target_average")

arm_summary <- summarize_groups(sample_slopes, c("target_label", "profiler", "assembly_arm"),
                                "response_slope", "descriptive_arm_slopes")

write.table(sample_slopes, file.path(outdir, "sample_arm_response_slopes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(paired_differences, file.path(outdir, "sample_paired_slope_differences.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(primary, file.path(outdir, "primary_assembly_effects.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(profiler_did, file.path(outdir, "secondary_profiler_differences.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pooled_summary, file.path(outdir, "secondary_pooled_assembly_effects.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(arm_summary, file.path(outdir, "descriptive_arm_slope_summaries.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

diagnostics <- data.frame(
  metric = c("status", "rows", "biological_samples", "targets", "profilers", "arms",
             "doses_per_cell", "sample_arm_slopes", "paired_differences",
             "bootstrap_replicates", "signflip_requested_replicates", "random_seed"),
  value = c("PASS", nrow(dat), length(unique(dat$sample_id)), length(unique(dat$target_label)),
            length(unique(dat$profiler)), length(unique(dat$assembly_arm)),
            length(expected_doses), nrow(sample_slopes), nrow(paired_differences),
            bootstrap_replicates, signflip_replicates, seed))
write.table(diagnostics, file.path(outdir, "sample_level_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "session_info.txt"))
seal_files <- c(normalizePath(input), list.files(outdir, full.names = TRUE))
seal_files <- seal_files[basename(seal_files) != "sample_level_analysis.sha256"]
status <- system2("sha256sum", seal_files,
                  stdout = file.path(outdir, "sample_level_analysis.sha256"))
if (!identical(status, 0L)) stop("Could not create SHA-256 seal.")
writeLines(c(paste("cohort", cohort_arg, sep = "\t"),
             paste("population", population_arg, sep = "\t"),
             paste("biological_samples", length(unique(dat$sample_id)), sep = "\t"),
             "inference_unit\tbiological_sample", "status\tPASS"),
           file.path(outdir, "SUCCESS"))
message("[PASS] Sample-level assembly sensitivity completed: ", outdir)
