#!/usr/bin/env Rscript
set.seed(41)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("assembly_sample_level."); dir.create(root, recursive = TRUE)
input <- file.path(root, "endpoints.tsv"); outdir <- file.path(root, "model")
grid <- expand.grid(sample_id = paste0("S", 1:8), target_label = c("Pana", "Pint"),
                    profiler = c("kraken2_bracken", "metaphlan4"),
                    assembly_arm = c("original", "clean"),
                    spike_fraction_target = c(.0001, .0005, .001, .005, .01, .05),
                    stringsAsFactors = FALSE)
grid$cohort <- "yachida"; grid$analysis_population <- "independent"
grid$condition <- ifelse(as.integer(sub("S", "", grid$sample_id)) <= 4, "Control", "CRC")
# Achieved fractions vary slightly by biological sample in real spike designs.
grid$spike_fraction_target <- grid$spike_fraction_target *
  (1 + as.integer(sub("S", "", grid$sample_id)) * 1e-5)
base <- ifelse(grid$profiler == "kraken2_bracken", .9, .7)
gain <- ifelse(grid$assembly_arm == "clean",
               ifelse(grid$target_label == "Pana", .12, .06), 0)
sample_effect <- (as.integer(sub("S", "", grid$sample_id)) - 4.5) / 200
grid$recovered_spike_signal <- (base + gain + sample_effect) * grid$spike_fraction_target +
  rnorm(nrow(grid), sd = 0.000002)
write.table(grid, input, sep = "\t", quote = FALSE, row.names = FALSE)
script <- file.path(repo, "analysis_v2/scripts/fit_assembly_sensitivity_sample_level.R")
result <- system2("Rscript", c(script, "--input", input, "--outdir", outdir,
                               "--bootstrap-replicates", "1000",
                               "--signflip-replicates", "1000"),
                  stdout = TRUE, stderr = TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))
stopifnot(file.exists(file.path(outdir, "SUCCESS")))
slopes <- read.delim(file.path(outdir, "sample_arm_response_slopes.tsv"))
diffs <- read.delim(file.path(outdir, "sample_paired_slope_differences.tsv"))
primary <- read.delim(file.path(outdir, "primary_assembly_effects.tsv"))
stopifnot(nrow(slopes) == 8 * 2 * 2 * 2, nrow(diffs) == 8 * 2 * 2, nrow(primary) == 4)
stopifnot(all(primary$n_samples == 8), all(primary$inference_method == "exact_sign_flip"))
stopifnot(all(is.finite(primary$mean)), all(primary$p_value > 0 & primary$p_value <= 1))
stopifnot(all(c("lower_95_bootstrap", "upper_95_bootstrap", "q_value_bh") %in% names(primary)))
cat("[PASS] sample-level assembly-sensitivity fixture\n")
