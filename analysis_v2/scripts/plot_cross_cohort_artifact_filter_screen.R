#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(scales)})
a<-commandArgs(TRUE);v<-function(x){i<-match(x,a);if(is.na(i)||i==length(a))stop("Missing ",x);a[i+1]}
input<-v("--input");out<-v("--outdir");d<-read.delim(input,check.names=FALSE)
stopifnot(all(c("train_cohort","test_cohort","analysis_population","profiler","score_threshold","off_target_reduction","target_recall_before","target_recall_after")%in%names(d)))
d$Transfer<-paste(tools::toTitleCase(d$train_cohort),"to",tools::toTitleCase(d$test_cohort));d$Population<-tools::toTitleCase(d$analysis_population);d$Profiler<-c(kraken2_bracken="Kraken2 + Bracken",metaphlan4="MetaPhlAn 4")[d$profiler]
d$recall_retained<-ifelse(d$target_recall_before>0,d$target_recall_after/d$target_recall_before,NA);d$threshold_label<-format(d$score_threshold,trim=TRUE)
write.table(d,file.path(out,"pareto_figure_source.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
p<-ggplot(d,aes(off_target_reduction,recall_retained,color=Transfer,group=Transfer))+
 geom_path(linewidth=.7)+geom_point(size=2)+geom_text(aes(label=threshold_label),size=2.4,nudge_y=.025,check_overlap=TRUE)+
 facet_grid(Population~Profiler)+geom_hline(yintercept=.95,linetype=2,color="grey40")+
 scale_x_continuous(labels=label_percent(),limits=c(0,1))+scale_y_continuous(labels=label_percent(),limits=c(0,1.08))+
 labs(x="Held-out off-target calls removed",y="Held-out target recall retained",color="Training to test",
 title="No general static blacklist transfers without target-recall loss",
 subtitle="Labels are recurrence thresholds; targets are not protected by an oracle",
 caption="Development-only post-fit screen. A valid filter must be learned pre-fit and followed by complete model/BH refitting.")+
 theme_bw(10)+theme(panel.grid.minor=element_blank(),legend.position="bottom",strip.background=element_rect(fill="grey94"),strip.text=element_text(face="bold"))
ggsave(file.path(out,"cross_cohort_filter_pareto.pdf"),p,width=10,height=6,device=cairo_pdf)
ggsave(file.path(out,"cross_cohort_filter_pareto.png"),p,width=10,height=6,dpi=320,bg="white")
writeLines(c("figure\tcross_cohort_filter_pareto","status\tPASS"),
           file.path(out,"FIGURE_SUCCESS"))
sealed<-c(normalizePath(input),file.path(out,"pareto_figure_source.tsv"),
          file.path(out,"cross_cohort_filter_pareto.pdf"),
          file.path(out,"cross_cohort_filter_pareto.png"),
          file.path(out,"FIGURE_SUCCESS"))
status<-system2("sha256sum",sealed,
                stdout=file.path(out,"artifact_filter_figure.sha256"))
if(!identical(status,0L))stop("Could not seal cross-cohort filter figure")
message("[PASS] cross-cohort filter screen figure")
