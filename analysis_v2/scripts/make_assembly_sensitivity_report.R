#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
value <- function(flag) { i <- match(flag,args); if(is.na(i)||i==length(args)) stop("Missing ",flag); args[[i+1L]] }
comparison <- value("--comparison-dir"); quantitative <- value("--quantitative-dir")
artificial <- value("--artificial-report"); outdir <- value("--outdir"); status <- value("--report-status")
if (!status %in% c("DEVELOPMENT_ONLY","DEFINITIVE")) stop("Invalid report status.")
required <- c(file.path(comparison,"SUCCESS"),file.path(comparison,"detection_arm_comparison.tsv"),
  file.path(comparison,"biomarker_arm_comparison.tsv"),file.path(quantitative,"SUCCESS"),
  file.path(quantitative,"primary_assembly_effects.tsv"),file.path(artificial,"SUCCESS"),
  file.path(artificial,"tables","off_target_call_ledger.tsv"))
if(any(!file.exists(required))||any(file.info(required)$size<=0)) stop("Assembly report input is incomplete.")
manifest <- read.delim(file.path(artificial,"provenance","report_manifest.tsv"),stringsAsFactors=FALSE)
source_status <- manifest$value[manifest$field=="status"]
if(status=="DEFINITIVE" && (!length(source_status)||source_status!="DEFINITIVE")) stop("Development evidence cannot produce a definitive report.")
if(dir.exists(outdir)&&length(list.files(outdir))) stop("OUTDIR must be new or empty.")
for(x in c("tables","figure_source","figures","provenance")) dir.create(file.path(outdir,x),recursive=TRUE,showWarnings=FALSE)
detection <- read.delim(file.path(comparison,"detection_arm_comparison.tsv"),check.names=FALSE)
biomarker <- read.delim(file.path(comparison,"biomarker_arm_comparison.tsv"),check.names=FALSE)
slopes <- read.delim(file.path(quantitative,"primary_assembly_effects.tsv"),check.names=FALSE)
if(!nrow(detection)||!nrow(biomarker)||!nrow(slopes)) stop("Empty assembly report table.")
detection$dose_percent <- 100*as.numeric(detection$spike_fraction_target)
biomarker$dose_percent <- 100*as.numeric(biomarker$spike_fraction_target)
biomarker$q_threshold <- as.numeric(biomarker$q_threshold)
primary_biomarker <- biomarker[abs(biomarker$q_threshold-.05)<1e-12,,drop=FALSE]
labels <- c(kraken2_bracken="Kraken2 + Bracken",metaphlan4="MetaPhlAn 4")
for(x in c("detection","primary_biomarker","slopes")) {
  frame <- get(x); frame$profiler_display <- unname(labels[frame$profiler]); assign(x,frame)
}
write.table(detection,file.path(outdir,"figure_source","detection_arm_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(primary_biomarker,file.path(outdir,"figure_source","biomarker_arm_comparison_q005.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(slopes,file.path(outdir,"figure_source","quantitative_slope_arm_comparison.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
file.copy(file.path(artificial,"tables","off_target_call_ledger.tsv"),file.path(outdir,"tables","off_target_call_ledger.tsv"),overwrite=TRUE)
file.copy(file.path(artificial,"tables","recurrent_off_target_taxa.tsv"),file.path(outdir,"tables","recurrent_off_target_taxa.tsv"),overwrite=TRUE)
theme_report <- theme_bw(base_size=10)+theme(legend.position="bottom",panel.grid.minor=element_blank())
p1 <- ggplot(detection,aes(dose_percent,as.numeric(clean_minus_original_detection_rate),color=profiler_display,group=profiler_display))+
  geom_hline(yintercept=0,linetype=2,color="grey50")+geom_line()+geom_point()+facet_wrap(~target_label)+
  labs(x="Implanted target fraction (%)",y="Clean minus original detection rate",color="Profiler")+theme_report
p2 <- ggplot(primary_biomarker,aes(dose_percent,as.numeric(clean_minus_original_target_effect),color=profiler_display,group=profiler_display))+
  geom_hline(yintercept=0,linetype=2,color="grey50")+geom_line()+geom_point()+facet_wrap(~target_label)+
  labs(x="Implanted target fraction (%)",y="Clean minus original target log2 effect",color="Profiler")+theme_report
p3 <- ggplot(primary_biomarker,aes(dose_percent,as.numeric(clean_minus_original_off_target_enriched_calls),color=profiler_display,group=profiler_display))+
  geom_hline(yintercept=0,linetype=2,color="grey50")+geom_line()+geom_point()+facet_wrap(~target_label)+
  labs(x="Implanted target fraction (%)",y="Clean minus original off-target calls",color="Profiler")+theme_report
p4 <- ggplot(slopes,aes(target_label,mean,color=profiler_display))+
  geom_hline(yintercept=0,linetype=2,color="grey50")+geom_pointrange(aes(ymin=lower_95_bootstrap,ymax=upper_95_bootstrap),position=position_dodge(width=.4))+
  labs(x="Target",y="Clean minus original response slope",color="Profiler")+theme_report
plots <- list(detection_change=p1,target_biomarker_effect_change=p2,off_target_burden_change=p3,quantitative_slope_change=p4)
for(name in names(plots)) { ggsave(file.path(outdir,"figures",paste0(name,".pdf")),plots[[name]],width=7,height=4.5); ggsave(file.path(outdir,"figures",paste0(name,".png")),plots[[name]],width=7,height=4.5,dpi=300) }
writeLines(c("# Draft assembly-sensitivity figure captions","",
  "Clean-minus-original differences compare alternative source assemblies for Pana and Pint in the same samples and doses. They are assembly-choice sensitivity estimates, not causal contamination effects.","",
  "Detection differences use paired binary outcomes; quantitative intervals use biological-sample bootstrap inference; biomarker and off-target differences are descriptive comparisons at BH q <= 0.05."),file.path(outdir,"captions.md"))
run_manifest <- data.frame(field=c("status","analysis","comparison","quantitative","artificial_report","created_at"),value=c(status,"assembly_choice_sensitivity",normalizePath(comparison),normalizePath(quantitative),normalizePath(artificial),format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")))
write.table(run_manifest,file.path(outdir,"provenance","report_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
files <- c(required,list.files(outdir,recursive=TRUE,full.names=TRUE)); rc <- system2("sha256sum",files,stdout=file.path(outdir,"provenance","report.sha256")); if(rc!=0) stop("Could not seal report.")
writeLines(c("report\tassembly_choice_sensitivity",paste0("report_status\t",status),"status\tPASS"),file.path(outdir,"SUCCESS"))
message("[PASS] Assembly-sensitivity report completed: ",outdir)
