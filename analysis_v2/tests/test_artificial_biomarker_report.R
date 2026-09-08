#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("artificial_report."); run <- file.path(root, "run"); out <- file.path(root, "report")
dir.create(file.path(run, "models"), recursive = TRUE); dir.create(file.path(run, "evaluation"))
writeLines("status\tPASS", file.path(run, "SUCCESS")); writeLines("status\tPASS", file.path(run, "models", "SUCCESS"))
writeLines("status\tPASS", file.path(run, "evaluation", "SUCCESS")); writeLines("status\tDEVELOPMENT_ONLY", file.path(run, "DEVELOPMENT_ONLY.txt"))
metrics <- expand.grid(target_label=c("Pana", "Pint"), assembly_arm=c("original", "clean"),
  profiler=c("kraken2_bracken", "metaphlan4"), condition=c("Control", "CRC"),
  spike_fraction_target=c(.0001, .001, .01), q_threshold=c(.05, .1), stringsAsFactors=FALSE)
metrics$cohort <- "yachida"; metrics$study <- "Study"; metrics$analysis_population <- "independent"
metrics$contrast <- paste0("spiked_vs_matched_baseline__background_", metrics$condition)
metrics$target_alias <- ifelse(metrics$target_label == "Pana", "Peptostreptococcus anaerobius", "Prevotella intermedia")
metrics$target_called <- as.integer(metrics$spike_fraction_target >= .001)
metrics$enriched_calls <- metrics$target_called + 1L; metrics$off_target_enriched_calls <- 1L
metrics$precision <- ifelse(metrics$target_called == 1, .5, 0); metrics$recall <- metrics$target_called
metrics$f1 <- ifelse(metrics$target_called == 1, 2/3, 0); metrics$target_effect <- log2(1 + 100 * metrics$spike_fraction_target)
metrics$target_q_value <- ifelse(metrics$target_called == 1, .01, .5)
metrics$target_effect_change_from_baseline <- metrics$target_effect
metrics$baseline_reference_kind <- "structural_null"; metrics$biomarker_set_jaccard_vs_baseline <- NA
write.table(metrics, file.path(run, "evaluation", "biomarker_propagation_metrics.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
calls <- data.frame(feature="target", effect=1, p_value=.01, q_value=.02, include=1)
write.table(calls, file.path(run, "models", "paired_da_results.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
script <- file.path(repo, "analysis_v2/scripts/make_artificial_biomarker_report.R")
result <- system2("Rscript", c(script, "--paired-run", run, "--outdir", out, "--report-status", "DEVELOPMENT_ONLY"), stdout=TRUE, stderr=TRUE)
status <- attr(result, "status"); if (is.null(status)) status <- 0L
if (status != 0L) stop(paste(result, collapse="\n"))
minimum <- read.delim(file.path(out, "tables", "minimum_detectable_dose.tsv"))
stopifnot(file.exists(file.path(out, "SUCCESS")), nrow(minimum) == 32,
          all(abs(minimum$minimum_detected_fraction - .001) < 1e-12),
          file.exists(file.path(out, "figures", "artificial_target_recall.pdf")),
          file.exists(file.path(out, "figures", "off_target_discovery_burden.png")),
          file.exists(file.path(out, "provenance", "report.sha256")))
cat("[PASS] artificial-biomarker report fixture\n")
