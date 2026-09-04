#!/usr/bin/env Rscript
set.seed(29)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("assembly_sensitivity_model.")
dir.create(root, recursive = TRUE)
input <- file.path(root, "endpoints.tsv"); outdir <- file.path(root, "model")
grid <- expand.grid(
  sample_id = paste0("S", 1:8), target_label = c("Pana", "Pint"),
  profiler = c("kraken2_bracken", "metaphlan4"),
  assembly_arm = c("original", "clean"),
  spike_fraction_target = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05),
  stringsAsFactors = FALSE)
grid$cohort <- "yachida"; grid$analysis_population <- "independent"
grid$condition <- ifelse(as.integer(sub("S", "", grid$sample_id)) <= 4L,
                         "Control", "CRC")
base_slope <- ifelse(grid$profiler == "kraken2_bracken", 0.9, 0.7)
clean_gain <- ifelse(grid$assembly_arm == "clean", 0.08, 0)
grid$recovered_spike_signal <- (base_slope + clean_gain) * grid$spike_fraction_target +
  rnorm(nrow(grid), sd = 0.00002)
write.table(grid, input, sep = "\t", quote = FALSE, row.names = FALSE)
script <- file.path(repo, "analysis_v2/scripts/fit_assembly_sensitivity.R")
result <- system2("Rscript", c(script, "--input", input, "--outdir", outdir),
                  stdout = TRUE, stderr = TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))
stopifnot(file.exists(file.path(outdir, "SUCCESS")))
slopes <- read.delim(file.path(outdir, "assembly_profiler_response_slopes.tsv"))
contrasts <- read.delim(file.path(outdir, "assembly_slope_contrasts.tsv"))
stopifnot(nrow(slopes) == 4L, nrow(contrasts) == 3L)
stopifnot(all(c("kraken2_bracken_clean_minus_original",
                "metaphlan4_clean_minus_original",
                "profiler_difference_in_differences") %in% contrasts$contrast))
cat("[PASS] assembly-sensitivity model synthetic fixture\n")
