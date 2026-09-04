#!/usr/bin/env Rscript
if (!requireNamespace("mgcv", quietly = TRUE)) {
  message("[SKIP] mgcv unavailable; run this test in the frozen analysis image")
  quit(status = 0)
}
args <- commandArgs(trailingOnly = FALSE)
this_file <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][[1]]))
repo <- normalizePath(file.path(dirname(this_file), "..", ".."))
script <- file.path(repo, "analysis_v2", "scripts", "fit_continuous_dose_response.R")
tmp <- tempfile("continuous_model."); dir.create(tmp)
input <- file.path(tmp, "endpoints.tsv"); outdir <- file.path(tmp, "model")
set.seed(73)
d <- expand.grid(sample_id = paste0("S", seq_len(24)),
  target_label = c("Fnuc", "Pmic", "Bfrag", "Pana"),
  spike_fraction_target = c(0.0001, 0.0005, 0.001, 0.005, 0.01),
  profiler = c("kraken2_bracken", "metaphlan4"), stringsAsFactors = FALSE)
d$cohort <- "synthetic"; d$condition <- rep(c("Control", "Adenoma", "CRC"), length.out = nrow(d))
d$analysis_population <- "independent"; d$assembly_arm <- "original"
sample_effect <- setNames(rnorm(24, 0, 0.00005), paste0("S", seq_len(24)))
target_effect <- setNames(rnorm(4, 0, 0.00003), c("Fnuc", "Pmic", "Bfrag", "Pana"))
slope <- ifelse(d$profiler == "kraken2_bracken", 0.92, 0.78)
d$recovered_spike_signal <- slope * d$spike_fraction_target + sample_effect[d$sample_id] +
  target_effect[d$target_label] + rnorm(nrow(d), 0, 0.00004 + 0.01 * d$spike_fraction_target)
write.table(d, input, sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(script, "--input", input, "--outdir", outdir,
  "--cohort", "synthetic", "--population", "independent"))
stopifnot(status == 0)
needed <- c("SUCCESS", "continuous_model_diagnostics.tsv", "profiler_response_slopes.tsv",
  "primary_profiler_slope_contrast.tsv", "categorical_dose_nonlinearity_test.tsv",
  "continuous_fitted_residuals.tsv", "continuous_model.sha256", "session_info.txt")
stopifnot(all(file.exists(file.path(outdir, needed))))
slopes <- read.delim(file.path(outdir, "profiler_response_slopes.tsv"))
stopifnot(nrow(slopes) == 2L, all(is.finite(slopes$response_slope)))
stopifnot(abs(slopes$response_slope[slopes$profiler == "kraken2_bracken"] - 0.92) < 0.08)
stopifnot(abs(slopes$response_slope[slopes$profiler == "metaphlan4"] - 0.78) < 0.08)
message("[PASS] continuous-model synthetic fixture")
