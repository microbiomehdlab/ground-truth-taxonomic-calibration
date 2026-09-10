#!/usr/bin/env Rscript
root <- tempfile("three_cohort_report."); dir.create(root)
cohorts <- c("yachida", "feng", "zeller")
artificial <- disease <- linkage <- character()
for (cohort in cohorts) {
  report <- file.path(root, cohort)
  for (directory in c("artificial/tables", "disease/tables", "linkage/figure_source"))
    dir.create(file.path(report, directory), recursive = TRUE)
  a <- data.frame(analysis_scope="pooled_primary",cohort=cohort,analysis_population="independent",
    profiler="kraken2_bracken",assembly_arm="original",q_threshold=.05,dose_percent_nominal=1,
    target_recall=.8,mean_precision=.7,context_sum_off_target_calls=2)
  d <- data.frame(cohort=cohort,analysis_population="independent",profiler="kraken2_bracken",
    contrast="CRC_vs_Control",q_threshold=.05,spike_fraction_target=.01,
    overall_baseline_retention=.85,median_jaccard=.8)
  l <- data.frame(analysis_scope="pooled_primary",cohort=cohort,analysis_population="independent",
    profiler="kraken2_bracken",target_label="Bfrag",q_threshold=.05,dose_percent_nominal=1,
    median_response_ratio=1.05,off_target_enriched_calls=2)
  write.table(a,file.path(report,"artificial/tables/artificial_biomarker_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(d,file.path(report,"disease/tables/disease_biomarker_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(l,file.path(report,"linkage/figure_source/calibration_biomarker_linkage.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  artificial <- c(artificial,file.path(report,"artificial")); disease <- c(disease,file.path(report,"disease")); linkage <- c(linkage,file.path(report,"linkage"))
}
out <- file.path(root,"out")
script <- file.path("analysis_v2","scripts","make_three_cohort_publication_report.R")
status <- system2("Rscript",c(script,"--artificial-reports",paste(artificial,collapse=","),
  "--disease-reports",paste(disease,collapse=","),"--linkage-reports",paste(linkage,collapse=","),
  "--report-status","DEVELOPMENT_ONLY","--outdir",out))
stopifnot(status==0,file.exists(file.path(out,"SUCCESS")),file.exists(file.path(out,"DEVELOPMENT_ONLY.txt")),
  file.exists(file.path(out,"figures","artificial_target_recall.pdf")),
  file.exists(file.path(out,"tables","three_cohort_disease_biomarkers.tsv")),
  file.exists(file.path(out,"provenance","publication_report.sha256")))
message("[PASS] three-cohort publication-report fixture")
