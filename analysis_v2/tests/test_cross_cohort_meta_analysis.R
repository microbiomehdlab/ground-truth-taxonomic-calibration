#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("meta_analysis."); dir.create(root, recursive = TRUE)
cohorts <- c("yachida", "feng", "zeller"); paths <- character()
for (i in seq_along(cohorts)) {
  rows <- expand.grid(dose_level = c("baseline", "dose_01"),
                      contrast = c("CRC_vs_Control", "Adenoma_vs_Control"),
                      feature = c("Target species", "Shared marker"), stringsAsFactors = FALSE)
  rows$cohort <- cohorts[i]; rows$study <- paste0("Study", i); rows$analysis_population <- "independent"
  rows$target_label <- "Pana"; rows$assembly_arm <- "original"; rows$profiler <- "metaphlan4"
  rows$spike_fraction_target <- ifelse(rows$dose_level == "baseline", 0, .001 * (1 + i * 1e-4))
  rows$effect <- ifelse(rows$feature == "Target species", .4 + i * .1, -.2 + i * .02)
  rows$standard_error <- .1; rows$p_value <- .01; rows$q_value <- .02
  rows$model_spec <- "primary_age_sex"; rows$include <- 1; rows$exclusion_reason <- ""
  path <- file.path(root, paste0(cohorts[i], ".tsv")); paths <- c(paths, path)
  write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
outdir <- file.path(root, "model")
script <- file.path(repo, "analysis_v2/scripts/meta_analyze_disease_biomarkers.R")
result <- system2("Rscript", c(script, "--inputs", paste(paths, collapse = ","),
  "--expected-cohorts", paste(cohorts, collapse = ","), "--outdir", outdir), stdout = TRUE, stderr = TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse = "\n"))
pooled <- read.delim(file.path(outdir, "random_effects_meta_analysis.tsv"))
stopifnot(file.exists(file.path(outdir, "SUCCESS")), nrow(pooled) == 8,
          all(pooled$cohort_count == 3), all(pooled$q_value >= 0 & pooled$q_value <= 1),
          all(is.finite(pooled$I2_pct)), all(pooled$direction_consistent == 1))
loo <- read.delim(file.path(outdir, "leave_one_cohort_out.tsv"))
stopifnot(nrow(loo) == 24)
cat("[PASS] cross-cohort meta-analysis fixture\n")
