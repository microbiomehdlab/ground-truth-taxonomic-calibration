#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=FALSE); f<-grep("^--file=",args,value=TRUE)
repo<-normalizePath(file.path(dirname(sub("^--file=","",f)),"../..")); root<-tempfile("pooled.");dir.create(root)
samples<-paste0("S",1:6); levels<-c("baseline","dose_01")
m<-expand.grid(sample_id=samples,dose_level=levels,stringsAsFactors=FALSE)
m$condition<-rep(c("Control","CRC"),each=3)[match(m$sample_id,samples)]
m$cohort<-"yachida";m$study<-"Study";m$analysis_population<-"independent";m$target_label<-"Pana";m$assembly_arm<-"original";m$profiler<-"kraken2_bracken"
m$profile_id<-paste(m$sample_id,m$dose_level);m$baseline_profile_id<-m$sample_id;m$spike_fraction_target<-ifelse(m$dose_level=="baseline",0,.01)
m$source_profile<-file.path(root,paste0(m$profile_id,".tsv"));m$target_taxon<-"Peptostreptococcus anaerobius";m$target_feature<-m$target_taxon;m$include<-1;m$exclusion_reason<-""
mp<-file.path(root,"manifest.tsv");write.table(m,mp,sep="\t",quote=FALSE,row.names=FALSE)
a<-do.call(rbind,lapply(seq_len(nrow(m)),function(i)data.frame(profiler=m$profiler[i],source_profile=m$source_profile[i],feature=m$target_taxon[i],abundance_fraction=.001+2*m$spike_fraction_target[i]+as.integer(sub("S","",m$sample_id[i]))*1e-6)))
ap<-file.path(root,"abundance.tsv");write.table(a,ap,sep="\t",quote=FALSE,row.names=FALSE);out<-file.path(root,"out")
r<-system2("Rscript",c(file.path(repo,"analysis_v2/scripts/fit_paired_biomarker_models.R"),"--profile-manifest",mp,"--abundance-long",ap,"--condition-mode","pooled","--outdir",out),stdout=TRUE,stderr=TRUE)
s<-attr(r,"status");if(is.null(s))s<-0L;if(s!=0L)stop(paste(r,collapse="\n"))
x<-read.delim(file.path(out,"paired_da_results.tsv"));stopifnot(file.exists(file.path(out,"SUCCESS")),length(unique(x$contrast))==1,unique(x$contrast)=="spiked_vs_matched_baseline__pooled_conditions",nrow(x[x$spike_fraction_target>0,])==1)
cat("[PASS] pooled paired biomarker-model fixture\n")
