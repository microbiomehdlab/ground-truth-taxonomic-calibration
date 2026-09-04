#!/usr/bin/env Rscript
if (!requireNamespace("mgcv", quietly = TRUE)) {
  message("[SKIP] mgcv unavailable; run this test in the frozen analysis image")
  quit(status = 0)
}

args <- commandArgs(trailingOnly = FALSE)
this_arg <- args[grepl("^--file=", args)]
this_file <- normalizePath(sub("^--file=", "", this_arg[[1]]))
repo <- normalizePath(file.path(dirname(this_file), "..", ".."))
script <- file.path(repo, "analysis_v2", "scripts", "fit_detection_dose_response.R")
tmp <- tempfile("detection_model.")
dir.create(tmp)
input <- file.path(tmp, "endpoints.tsv")
outdir <- file.path(tmp, "model")

set.seed(41)
samples <- paste0("S", seq_len(18))
targets <- c("Fnuc", "Pmic", "Bfrag")
doses <- c(0, 0.0001, 0.001, 0.01)
profilers <- c("kraken2_bracken", "metaphlan4")
d <- expand.grid(sample_id = samples, target_label = targets,
                 spike_fraction_target = doses, profiler = profilers,
                 stringsAsFactors = FALSE)
d$cohort <- "synthetic"
d$study <- "fixture"
d$condition <- rep(c("Control", "CRC"), length.out = nrow(d))
d$analysis_population <- "independent"
d$assembly_arm <- "original"
eta <- -2 + 1.4 * (d$spike_fraction_target > 0) +
  1.2 * (d$spike_fraction_target >= 0.001) + 0.35 * (d$profiler == "metaphlan4")
d$observed_detected_native_nonzero <- rbinom(nrow(d), 1, plogis(eta))
write.table(d, input, sep = "\t", quote = FALSE, row.names = FALSE)

status <- system2("Rscript", c(script, "--input", input, "--outdir", outdir,
                                "--cohort", "synthetic", "--population", "independent"))
stopifnot(status == 0)
needed <- c("SUCCESS", "detection_model_diagnostics.tsv",
            "detection_probability_curve.tsv", "profiler_contrasts_by_dose.tsv",
            "profiler_by_dose_primary_test.txt", "profiler_by_dose_primary_test.tsv",
            "detection_model.sha256", "session_info.txt")
stopifnot(all(file.exists(file.path(outdir, needed))))
curve <- read.delim(file.path(outdir, "detection_probability_curve.tsv"))
stopifnot(nrow(curve) == length(doses) * length(profilers))
message("[PASS] detection-model synthetic fixture")
