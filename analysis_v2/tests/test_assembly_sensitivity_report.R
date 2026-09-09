#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=FALSE); file_arg <- grep("^--file=",args,value=TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=","",file_arg)),"../.."))
root <- tempfile("assembly_report."); comparison <- file.path(root,"comparison"); quantitative <- file.path(root,"quantitative"); artificial <- file.path(root,"artificial"); out <- file.path(root,"out")
for(x in c(comparison,quantitative,file.path(artificial,"tables"),file.path(artificial,"provenance"))) dir.create(x,recursive=TRUE,showWarnings=FALSE)
writeLines("status\tPASS",file.path(comparison,"SUCCESS")); writeLines("status\tPASS",file.path(quantitative,"SUCCESS")); writeLines("status\tPASS",file.path(artificial,"SUCCESS"))
detection <- expand.grid(target_label=c("Pana","Pint"),profiler=c("kraken2_bracken","metaphlan4"),spike_fraction_target=c(.001,.01),stringsAsFactors=FALSE)
detection$samples <- 10; detection$original_detection_rate <- .5; detection$clean_detection_rate <- .7; detection$clean_minus_original_detection_rate <- .2; detection$clean_only_detections <- 2; detection$original_only_detections <- 0; detection$exact_mcnemar_p_value <- .5
write.table(detection,file.path(comparison,"detection_arm_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
biomarker <- expand.grid(target_label=c("Pana","Pint"),profiler=c("kraken2_bracken","metaphlan4"),contrast="pooled",spike_fraction_target=c(.001,.01),q_threshold=c(.05,.1),stringsAsFactors=FALSE)
for(field in c("target_called","target_effect","off_target_enriched_calls","precision")) { biomarker[[paste0("original_",field)]] <- 1; biomarker[[paste0("clean_",field)]] <- 2; biomarker[[paste0("clean_minus_original_",field)]] <- 1 }
write.table(biomarker,file.path(comparison,"biomarker_arm_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
slopes <- expand.grid(target_label=c("Pana","Pint"),profiler=c("kraken2_bracken","metaphlan4"),stringsAsFactors=FALSE); slopes$mean <- .1; slopes$lower_95_bootstrap <- -.1; slopes$upper_95_bootstrap <- .3
write.table(slopes,file.path(quantitative,"primary_assembly_effects.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(feature="offtarget"),file.path(artificial,"tables","off_target_call_ledger.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(feature="offtarget"),file.path(artificial,"tables","recurrent_off_target_taxa.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(field="status",value="DEVELOPMENT_ONLY"),file.path(artificial,"provenance","report_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
script <- file.path(repo,"analysis_v2/scripts/make_assembly_sensitivity_report.R")
result <- system2("Rscript",c(script,"--comparison-dir",comparison,"--quantitative-dir",quantitative,"--artificial-report",artificial,"--outdir",out,"--report-status","DEVELOPMENT_ONLY"),stdout=TRUE,stderr=TRUE)
rc <- attr(result,"status"); if(is.null(rc)) rc <- 0L; if(rc!=0) stop(paste(result,collapse="\n"))
stopifnot(file.exists(file.path(out,"SUCCESS")),file.exists(file.path(out,"figures","detection_change.pdf")),file.exists(file.path(out,"figures","off_target_burden_change.png")),file.exists(file.path(out,"provenance","report.sha256")))
cat("[PASS] assembly-sensitivity report fixture\n")
