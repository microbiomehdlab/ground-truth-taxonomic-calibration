#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=FALSE); file_arg<-grep("^--file=",args,value=TRUE)
repo<-normalizePath(file.path(dirname(sub("^--file=","",file_arg)),"../.."))
root<-tempfile("linkage."); dir.create(root); out<-file.path(root,"out")
conditions<-c("Control","CRC"); targets<-c("Pana","Pint"); profilers<-c("kraken2_bracken","metaphlan4")
doses<-c(.0001,.001,.01); samples<-paste0("S",1:3)
e<-expand.grid(condition=conditions,target_label=targets,profiler=profilers,spike_fraction_target=doses,
 sample_id=samples,stringsAsFactors=FALSE)
e$cohort<-"yachida";e$study<-"Study";e$analysis_population<-"independent";e$assembly_arm<-"original"
e$response_ratio<-ifelse(e$profiler=="metaphlan4",2,.9);e$recovered_spike_signal<-e$response_ratio*e$spike_fraction_target
e$signed_reference_error<-(e$response_ratio-1)*e$spike_fraction_target;e$absolute_reference_error<-abs(e$signed_reference_error)
endpoints<-file.path(root,"endpoints.tsv");write.table(e,endpoints,sep="\t",quote=FALSE,row.names=FALSE)
m<-unique(e[c("cohort","study","analysis_population","condition","target_label","assembly_arm","profiler","spike_fraction_target")])
m$contrast<-paste0("spiked_vs_matched_baseline__background_",m$condition);m$q_threshold<-.05
m$target_called<-as.integer(m$spike_fraction_target>=.001);m$target_effect<-log2(1+100*m$spike_fraction_target)
m$target_q_value<-ifelse(m$target_called,.01,.5);m$precision<-m$target_called;m$off_target_enriched_calls<-ifelse(m$profiler=="kraken2_bracken",2,0)
metrics<-file.path(root,"metrics.tsv");write.table(m,metrics,sep="\t",quote=FALSE,row.names=FALSE)
script<-file.path(repo,"analysis_v2/scripts/link_calibration_to_biomarkers.R")
result<-system2("Rscript",c(script,"--endpoints",endpoints,"--biomarker-metrics",metrics,"--outdir",out,"--analysis-status","DEVELOPMENT_ONLY"),stdout=TRUE,stderr=TRUE)
status<-attr(result,"status");if(is.null(status))status<-0L;if(status!=0L)stop(paste(result,collapse="\n"))
linked<-read.delim(file.path(out,"figure_source","calibration_biomarker_linkage.tsv"))
stopifnot(nrow(linked)==nrow(m),all(linked$biological_samples==3),
 file.exists(file.path(out,"tables","calibration_biomarker_associations.tsv")),
 file.exists(file.path(out,"figures","calibration_vs_target_effect.pdf")),
 file.exists(file.path(out,"provenance","linkage.sha256")),file.exists(file.path(out,"SUCCESS")))
cat("[PASS] calibration-to-biomarker linkage fixture\n")
