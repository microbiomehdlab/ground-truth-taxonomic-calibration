#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(scales)})
a<-commandArgs(TRUE);v<-function(x){i<-match(x,a);if(is.na(i)||i==length(a))stop("Missing ",x);a[i+1]}
root<-v("--calibration-root");out<-v("--outdir")
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("OUTDIR must be new or empty")
ledger<-read.delim(file.path(root,"transfer_ledger.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
if(!nrow(ledger)||!all(c("transfer","directory","validation_cohort")%in%names(ledger)))stop("Invalid transfer ledger")
abundance<-list();biomarker<-list();coefficients<-list()
for(i in seq_len(nrow(ledger))){
 d<-ledger$directory[i]
 ap<-file.path(d,"abundance_restoration_metrics.tsv")
 bp<-file.path(d,"restoration","calibration_restoration_summary.tsv")
 cp<-file.path(d,"calibration_coefficients.tsv")
 if(!file.exists(ap)||!file.exists(bp)||!file.exists(cp))stop("Incomplete transfer: ",ledger$transfer[i])
 x<-read.delim(ap,check.names=FALSE);x$transfer<-ledger$transfer[i];abundance[[i]]<-x
 y<-read.delim(bp,check.names=FALSE);y$transfer<-ledger$transfer[i];biomarker[[i]]<-y
 z<-read.delim(cp,check.names=FALSE);z$transfer<-ledger$transfer[i];coefficients[[i]]<-z
}
ad<-do.call(rbind,abundance);bd<-do.call(rbind,biomarker);cd<-do.call(rbind,coefficients)
needed.a<-c("validation_cohort","analysis_population","profiler","relative_mae_restoration")
needed.b<-c("validation_cohort","analysis_population","profiler","spike_fraction_target",
 "baseline_biomarkers_rescued","baseline_biomarkers_harmed","induced_calls_removed","induced_calls_created",
 "relative_mean_effect_restoration")
if(!all(needed.a%in%names(ad))||!all(needed.b%in%names(bd)))stop("Restoration columns incomplete")
prof<-c(kraken2_bracken="Kraken2 + Bracken",metaphlan4="MetaPhlAn 4")
decorate<-function(d){d$Cohort<-tools::toTitleCase(d$validation_cohort);d$Population<-tools::toTitleCase(d$analysis_population);d$Profiler<-unname(prof[d$profiler]);d}
ad<-decorate(ad);bd<-decorate(bd);cd<-decorate(cd);ad<-ad[is.finite(ad$relative_mae_restoration),];bd<-bd[is.finite(bd$relative_mean_effect_restoration),]
dir.create(file.path(out,"figures"),recursive=TRUE);dir.create(file.path(out,"figure_source"),recursive=TRUE)
theme_set(theme_bw(11)+theme(strip.text=element_text(face="bold"),legend.position="bottom",plot.title.position="plot"))
active<-cd[cd$correction_active==1 & is.finite(cd$slope_fraction_per_implanted_fraction) & cd$slope_fraction_per_implanted_fraction>0,]
if(!nrow(active))stop("No active calibration coefficients available for plotting")
panel<-interaction(active$transfer,active$analysis_population,active$profiler,drop=TRUE)
keep<-unlist(lapply(split(seq_len(nrow(active)),panel),function(ii){
 ii[order(active$slope_fraction_per_implanted_fraction[ii],decreasing=TRUE)][seq_len(min(15,length(ii)))]
}),use.names=FALSE)
atlas<-active[sort(unique(keep)),]
atlas$Feature<-factor(atlas$feature,levels=rev(unique(atlas$feature[order(atlas$slope_fraction_per_implanted_fraction,decreasing=TRUE)])))
atlas$Perturbation<-atlas$perturbation
p0<-ggplot(atlas,aes(Perturbation,Feature,fill=slope_fraction_per_implanted_fraction))+
 geom_tile(color="white",linewidth=.2)+facet_grid(Population~interaction(Cohort,Profiler,sep="\n"),scales="free_y",space="free_y")+
 scale_fill_viridis_c(option="C",labels=label_number(accuracy=.01))+
 labs(title="Training cohorts reveal recurrent target-to-taxon abundance inflation",subtitle="Top active held-out calibration coefficients; community perturbations are estimated jointly",x="Spike perturbation",y="Non-target feature",fill="Inflation at\n1% spike (pp)")+
 theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="right")
p1<-ggplot(ad,aes(Population,relative_mae_restoration,fill=Population))+
 geom_hline(yintercept=0,linetype=2,color="grey45")+geom_boxplot(outlier.shape=NA,alpha=.65)+geom_jitter(width=.15,alpha=.18,size=.6)+
 facet_grid(Cohort~Profiler)+scale_y_continuous(labels=label_percent())+
 labs(title="Does training-cohort calibration restore held-out abundance profiles?",subtitle="Positive values move non-target abundances toward the ideal dilution reference",x=NULL,y="Relative abundance-error restoration",fill="Population")

outcome_rows<-function(field,label,sign=1){
 x<-bd[,c("Cohort","Population","Profiler","spike_fraction_target")]
 x$outcome<-label;x$count<-sign*bd[[field]];x
}
long<-rbind(
 outcome_rows("baseline_biomarkers_rescued","Baseline biomarkers rescued"),
 outcome_rows("baseline_biomarkers_harmed","Baseline biomarkers harmed",-1),
 outcome_rows("induced_calls_removed","Induced calls removed"),
 outcome_rows("induced_calls_created","Induced calls created",-1))
long$Dose<-100*long$spike_fraction_target
p2<-ggplot(long,aes(Dose,count,color=outcome,group=outcome))+
 geom_hline(yintercept=0,color="grey45")+stat_summary(fun=median,geom="line",linewidth=.7)+stat_summary(fun=median,geom="point",size=1.7)+
 facet_grid(Population~interaction(Cohort,Profiler,sep="\n"),scales="free_y")+
 labs(title="Held-out calibration can rescue or damage biomarker calls",subtitle="Positive values are desired outcomes; negative values are newly introduced errors",x="Implanted target fraction (%)",y="Median biomarker calls per context",color="Outcome")

p3<-ggplot(bd,aes(100*spike_fraction_target,relative_mean_effect_restoration,color=Profiler))+
 geom_hline(yintercept=0,linetype=2,color="grey45")+stat_summary(fun=median,geom="line")+stat_summary(fun=median,geom="point")+
 facet_grid(Population~Cohort)+scale_y_continuous(labels=label_percent())+
 labs(title="Calibration is judged against the original disease-effect estimates",subtitle="Positive values indicate movement toward the unspiked CRC-versus-Control effect",x="Implanted target fraction (%)",y="Relative disease-effect restoration",color="Profiler")
for(s in list(list(p=p0,n="distortion_coefficient_atlas",w=14,h=9),list(p=p1,n="abundance_error_restoration",w=11,h=7),list(p=p2,n="biomarker_call_restoration",w=14,h=8),list(p=p3,n="disease_effect_restoration",w=10,h=7))){
 ggsave(file.path(out,"figures",paste0(s$n,".pdf")),s$p,width=s$w,height=s$h,device=cairo_pdf)
 ggsave(file.path(out,"figures",paste0(s$n,".png")),s$p,width=s$w,height=s$h,dpi=300,bg="white")
}
write.table(ad,file.path(out,"figure_source","abundance_restoration.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(bd,file.path(out,"figure_source","biomarker_restoration.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(atlas,file.path(out,"figure_source","distortion_coefficient_atlas.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
writeLines(c("status=DEVELOPMENT_ONLY","automatic_feature_removal=NO","heldout_training_separation=YES"),file.path(out,"DEVELOPMENT_ONLY.txt"))
writeLines(c("# Calibration-restoration figures","","The distortion atlas shows the strongest positive abundance-inflation coefficients learned outside each displayed validation cohort.","All panels evaluate rules learned outside the displayed validation cohort.","Positive restoration moves the perturbed result toward its unspiked reference.","Direct implanted targets are protected from subtraction. Results remain development-only."),file.path(out,"FIGURE_GUIDE.md"))
files<-c(normalizePath(file.path(root,"transfer_ledger.tsv")),list.files(out,recursive=TRUE,full.names=TRUE));files<-files[!grepl("SUCCESS$|calibration_restoration_figures.sha256$",files)]
rc<-system2("sha256sum",files,stdout=file.path(out,"calibration_restoration_figures.sha256"));if(!identical(rc,0L))stop("Could not seal figures")
writeLines("status=PASS",file.path(out,"SUCCESS"));message("[PASS] Calibration-restoration figures: ",normalizePath(out))
