#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("report."); run <- file.path(root, "run"); out <- file.path(root, "report")
dir.create(file.path(run, "models"), recursive = TRUE); dir.create(file.path(run, "evaluation"))
writeLines("status\tPASS", file.path(run, "SUCCESS")); writeLines("status\tPASS", file.path(run, "models", "SUCCESS"))
writeLines("status\tPASS", file.path(run, "evaluation", "SUCCESS")); writeLines("status\tDEVELOPMENT_ONLY", file.path(run, "DEVELOPMENT_ONLY.txt"))
metrics <- expand.grid(profiler = c("kraken2_bracken", "metaphlan4"), contrast = c("CRC_vs_Control", "Adenoma_vs_Control"),
  spike_fraction_target = c(.001, .01), q_threshold = c(.05, .1), stringsAsFactors = FALSE)
metrics$cohort <- "yachida"; metrics$study <- "Study"; metrics$analysis_population <- "independent"
metrics$target_label <- "Pana"; metrics$assembly_arm <- "clean"; metrics$target_alias <- "Target species"
metrics$baseline_biomarkers <- 2; metrics$dose_biomarkers <- 2; metrics$retained_biomarkers <- 1
metrics$lost_biomarkers <- 1; metrics$gained_biomarkers <- 1; metrics$baseline_retention_rate <- .5
metrics$dose_overlap_fraction <- .5; metrics$direction_flips_among_retained <- 0
metrics$median_abs_effect_change_baseline_biomarkers <- .1; metrics$max_abs_effect_change_baseline_biomarkers <- .2
metrics$target_significant <- 1; metrics$target_effect <- .4
metrics$target_q_value <- .01; metrics$target_effect_change_from_baseline <- .2
metrics$baseline_reference_kind <- "observed_calls"; metrics$biomarker_set_jaccard_vs_baseline <- .5
write.table(metrics, file.path(run, "evaluation", "disease_biomarker_propagation_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
models <- data.frame(cohort="yachida", analysis_population="independent", target_label="Pana", assembly_arm="clean",
 profiler="metaphlan4", effect=.4, standard_error=.1, lower_95=.2, upper_95=.6, model_spec="primary_age_sex",
 n_samples=10, covariates_used="age", covariates_omitted="sex_invariant")
write.table(models, file.path(run, "models", "primary_disease_da_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
script <- file.path(repo, "analysis_v2/scripts/make_disease_biomarker_report.R")
result <- system2("Rscript", c(script, "--disease-run", run, "--outdir", out, "--report-status", "DEVELOPMENT_ONLY"), stdout=TRUE, stderr=TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse="\n"))
stopifnot(file.exists(file.path(out, "SUCCESS")), file.exists(file.path(out, "figures", "baseline_biomarker_retention.pdf")),
          file.exists(file.path(out, "figures", "biomarker_set_stability.png")),
          nrow(read.delim(file.path(out, "tables", "disease_biomarker_summary.tsv"))) == 16)
cat("[PASS] disease-biomarker report fixture\n")
