#!/usr/bin/env Rscript
if (!requireNamespace("ggplot2", quietly = TRUE)) { message("[SKIP] ggplot2 unavailable"); quit(status = 0) }
args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", args, value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
root <- tempfile("combined_calibration."); models <- file.path(root, "models"); out <- file.path(root, "out")
for (cohort in c("feng", "zeller")) for (population in c("independent", "community")) {
  d <- file.path(models, paste0("detection_", cohort, "_", population)); q <- file.path(models, paste0("continuous_", cohort, "_", population))
  dir.create(d, recursive = TRUE); dir.create(q, recursive = TRUE)
  writeLines("status\tPASS", file.path(d, "SUCCESS")); writeLines("status\tPASS", file.path(q, "SUCCESS"))
  curve <- expand.grid(profiler=c("kraken2_bracken","metaphlan4"),spike_fraction_target=c(0,.001,.01)); curve$predicted_detection_probability <- c(.1,.2,.5,.7,.8,.9); curve$lower_95_descriptive <- pmax(0,curve$predicted_detection_probability-.05); curve$upper_95_descriptive <- pmin(1,curve$predicted_detection_probability+.05)
  contrast <- data.frame(spike_fraction_target=c(0,.001),odds_ratio=c(1,2))
  slopes <- data.frame(profiler=c("kraken2_bracken","metaphlan4"),response_slope=c(.9,1.2),lower_95=c(.8,1.1),upper_95=c(1,1.3))
  slope_contrast <- data.frame(term="interaction",difference_response_slope_metaphlan_minus_bracken=.3)
  residual <- expand.grid(profiler=c("kraken2_bracken","metaphlan4"),spike_fraction_target=c(.001,.01)); residual$residual_median <- c(-.1,.1,-.05,.05)
  diag <- data.frame(metric="converged",value=1)
  write.table(curve,file.path(d,"detection_probability_curve.tsv"),sep="\t",quote=FALSE,row.names=FALSE); write.table(contrast,file.path(d,"profiler_contrasts_by_dose.tsv"),sep="\t",quote=FALSE,row.names=FALSE); write.table(diag,file.path(d,"detection_model_diagnostics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(slopes,file.path(q,"profiler_response_slopes.tsv"),sep="\t",quote=FALSE,row.names=FALSE); write.table(slope_contrast,file.path(q,"primary_profiler_slope_contrast.tsv"),sep="\t",quote=FALSE,row.names=FALSE); write.table(residual,file.path(q,"continuous_residuals_by_dose.tsv"),sep="\t",quote=FALSE,row.names=FALSE); write.table(diag,file.path(q,"continuous_model_diagnostics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}
script <- file.path(repo,"analysis_v2/scripts/make_combined_calibration_report.R")
result <- system2("Rscript",c(script,"--models-root",models,"--outdir",out,"--report-status","DEVELOPMENT_ONLY"),stdout=TRUE,stderr=TRUE); rc <- attr(result,"status"); if(is.null(rc)) rc <- 0L; if(rc!=0) stop(paste(result,collapse="\n"))
stopifnot(file.exists(file.path(out,"SUCCESS")),file.exists(file.path(out,"tables","detection_profiler_contrasts.tsv")),file.exists(file.path(out,"figures","detection_probability.pdf")),file.exists(file.path(out,"figures","quantitative_response_slopes.png")),file.exists(file.path(out,"provenance","report.sha256")))
cat("[PASS] combined calibration report fixture\n")
