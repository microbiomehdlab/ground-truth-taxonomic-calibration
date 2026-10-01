#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
value <- function(flag) {
  i <- match(flag,args)
  if (is.na(i) || i == length(args)) stop('Required: ',flag)
  args[[i+1L]]
}
input <- value('--summary-dir')
out <- value('--outdir')
if (dir.exists(out) || file.exists(out)) stop('Choose a fresh figure directory')
status <- jsonlite::fromJSON(file.path(input,'summary.json'))
if (status$status != 'PASS_SUMMARY_TABLES') stop('Summary stage failed')
data <- read.delim(file.path(input,'recovery_summary.tsv'),check.names=FALSE)
differences <- read.delim(file.path(input,'paired_method_differences.tsv'),check.names=FALSE)
if (!nrow(data) || any(!is.finite(data$recovery_ratio_median))) stop('Invalid summaries')
is_mpa <- identical(status$profiler, 'metaphlan4')
profiler_title <- if(is_mpa) 'MetaPhlAn' else 'Bracken'
method_order <- if(is_mpa) c('Read-fraction reference / sensitivity','Genome-equivalent reference / primary') else c('Native fraction / old reference','All input pairs / read reference')
data$method <- factor(data$method,levels=method_order)
taxa <- c('Bfrag','Csym','Dpne','Fnuc','Hhat','Pmic','Pana','Psto','Porp','Pint')
unknown <- setdiff(unique(data$target_label),taxa)
if(length(unknown)) stop('Unknown taxa: ',paste(unknown,collapse=','))
data$target_label <- factor(data$target_label,levels=taxa)
differences$target_label <- factor(differences$target_label,levels=taxa)
data$condition <- factor(data$condition,levels=c('Control','Adenoma','CRC'))
differences$condition <- factor(differences$condition,levels=c('Control','Adenoma','CRC'))
dir.create(out,recursive=TRUE)
colors <- setNames(c('#737373','#0072B2'),method_order)
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
  legend.position='bottom',strip.background=element_rect(fill='#F3F3F3'),
  plot.title=element_text(face='bold'),plot.subtitle=element_text(size=10)))
save_plot <- function(plot,name) {
  ggsave(file.path(out,paste0(name,'.pdf')),plot,width=22,height=9,limitsize=FALSE)
  ggsave(file.path(out,paste0(name,'.png')),plot,width=22,height=9,dpi=180,limitsize=FALSE)
}
for (cohort in unique(data$cohort)) for (population in unique(data$analysis_population)) {
  d <- data[data$cohort==cohort & data$analysis_population==population,,drop=FALSE]
  delta <- differences[differences$cohort==cohort & differences$analysis_population==population,,drop=FALSE]
  if (!nrow(d)) next
  xlabel <- if(population=='community') 'Nominal total community addition (%)' else 'Nominal individual addition (%)'
  context <- paste(if(isTRUE(status$fixture)) 'SYNTHETIC TEST DATA' else 'Observed data', cohort,if(population=='community') 'community spikes' else 'individual spikes',sep=' | ')
  subtitle <- 'Points: sample medians; bars: sample Q1-Q3. All observations retained. Vertical axis uses a signed pseudo-log scale.'
  base <- ggplot(d,aes(x=100*nominal_total_dose,color=method,group=method))+
    facet_grid(condition~target_label,drop=FALSE)+scale_x_log10()+
    scale_color_manual(values=colors,drop=FALSE)+labs(x=xlabel,color=NULL,
      subtitle=subtitle,caption='Source: recovery_summary.tsv. IQR describes sample variation; it is not a confidence interval.')
  recovery <- base+geom_hline(yintercept=1,color='#333333',linetype='dashed')+
    geom_linerange(aes(ymin=recovery_ratio_q1,ymax=recovery_ratio_q3),linewidth=.35)+
    geom_line(aes(y=recovery_ratio_median),linewidth=.4)+geom_point(aes(y=recovery_ratio_median),size=1.4)+
    scale_y_continuous(trans=scales::pseudo_log_trans(base=10,sigma=.2))+
    labs(title=paste(profiler_title,'recovery ratio:',context),y='Recovered / implanted signal (ideal = 1)')
  error <- base+geom_linerange(aes(ymin=absolute_relative_error_q1,ymax=absolute_relative_error_q3),linewidth=.35)+
    geom_line(aes(y=absolute_relative_error_median),linewidth=.4)+geom_point(aes(y=absolute_relative_error_median),size=1.4)+
    scale_y_continuous(trans=scales::pseudo_log_trans(base=10,sigma=.2))+
    labs(title=paste(profiler_title,'absolute relative error:',context),y='Absolute error / implanted signal (ideal = 0)')
  delta_plot <- ggplot(delta,aes(x=100*nominal_total_dose,y=paired_error_delta_median))+
    geom_hline(yintercept=0,linetype='dashed',color='#333333')+
    geom_linerange(aes(ymin=paired_error_delta_q1,ymax=paired_error_delta_q3),color='#0072B2',linewidth=.35)+
    geom_line(color='#0072B2',linewidth=.4)+geom_point(color='#0072B2',size=1.4)+
    facet_grid(condition~target_label,drop=FALSE)+scale_x_log10()+
    scale_y_continuous(trans=scales::pseudo_log_trans(base=10,sigma=.2))+
    labs(title=paste('Paired change in recovery error:',context),x=xlabel,
      y=if(is_mpa) 'Genome-equivalent error minus read-reference error' else 'All-input error minus old error',subtitle='Negative: lower error using the primary reference. Points: median paired difference; bars: Q1-Q3.',
      caption='Source: paired_method_differences.tsv. Descriptive comparison; no inferential significance claim.')
  prefix <- paste(cohort,population,sep='_')
  save_plot(recovery,paste0(prefix,'_recovery_ratio'))
  save_plot(error,paste0(prefix,'_absolute_relative_error'))
  save_plot(delta_plot,paste0(prefix,'_paired_error_change'))
}
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
writeLines(c('status\tPASS_PLOTS','intervals\tsample_IQR','y_axis\tpseudo_log_signed'),file.path(out,'SUCCESS'))
message('[PASS] ',profiler_title,' comparison figures: ',out)
