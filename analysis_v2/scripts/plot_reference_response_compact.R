# Pilot-only descriptive synthesis. No new statistical fits.
args <- commandArgs(TRUE); stopifnot(length(args)==4)
report <- normalizePath(args[1]); audit <- normalizePath(args[2])
out <- args[3]; style <- normalizePath(args[4])
if(file.exists(out)) stop('Fresh figure output required')
verify <- function(root) {
  old <- getwd();on.exit(setwd(old));setwd(root)
  if(system2('sha256sum',c('-c','--quiet','SHA256SUMS'))!=0) stop('Input checksum failure')
}
verify(report);verify(audit);source(style);library(ggplot2)
read <- function(root,name) read.delim(file.path(root,name),stringsAsFactors=FALSE,check.names=FALSE)
d <- read(report,'target_comparisons.tsv');ne <- read(audit,'input_summary.tsv')
key <- function(x) paste(x$context_id,x$target_label,sep='|')
stopifnot(nrow(d)==7200,length(unique(d$context_id))==720,!anyDuplicated(key(d)),
          !anyDuplicated(key(ne)),setequal(key(ne),key(d[d$observed_estimable==FALSE,])))
stopifnot(all(tolower(as.character(ne$saved_status_agrees))=='true'),
          all(ne$constant_kind[ne$profiler=='metaphlan4']=='ALL_ZERO'),
          all(ne$native_unique_values[ne$profiler=='kraken2_bracken']==2))
tools <- c(kraken2_bracken='Kraken2/Bracken',metaphlan4='MetaPhlAn 4')
x <- d[d$n==20,];stopifnot(nrow(x)==3600)
summary <- function(data,by) {
  cells <- split(data,interaction(data[,by,drop=FALSE],drop=TRUE))
  do.call(rbind,lapply(cells,function(cell) {
    common <- cell[1,by,drop=FALSE]
    do.call(rbind,lapply(c('observed','reference'),function(arm) {
      data.frame(common,arm=arm,count=sum(cell[[paste0(arm,'_positive_discovery')]]),
                 denominator=nrow(cell),frequency=mean(cell[[paste0(arm,'_positive_discovery')]]))
    }))
  }))
}
pooled <- summary(x,c('profiler','nominal_total_dose'))
low <- summary(x[x$nominal_total_dose==.0001,],c('profiler','target_label'))
stopifnot(all(pooled$denominator==600),all(low$denominator==60))
decorate <- function(z) {
  z$Profiler <- factor(tools[z$profiler],levels=MANUSCRIPT_PROFILERS)
  z$Arm <- factor(z$arm,levels=c('observed','reference'),labels=c('Observed','Reference'))
  z
}
pooled <- decorate(pooled);low <- decorate(low)
pooled$Dose <- factor(pooled$nominal_total_dose,levels=c(.0001,.001,.01),labels=c('0.001','0.01','0.1'))
reason <- data.frame(profiler=names(tools),denominator=600L,count=c(
  sum(ne$profiler=='kraken2_bracken' & ne$n==20 & ne$nominal_total_dose==.0001),
  sum(ne$profiler=='metaphlan4' & ne$n==20 & ne$nominal_total_dose==.0001)),
  explanation=c('Perfect fixed fits','All-zero inputs'))
reason$Profiler <- factor(tools[reason$profiler],levels=MANUSCRIPT_PROFILERS)
reason$frequency <- reason$count/reason$denominator
dir.create(out,recursive=TRUE);out <- normalizePath(out)
dir.create(file.path(out,'source_data'));dir.create(file.path(out,'provenance'))
write <- function(z,name) write.table(z,file.path(out,'source_data',name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write(pooled,'pooled_frequencies.tsv');write(low,'low_dose_target_frequencies.tsv');write(reason,'ne_reasons.tsv')
a <- ggplot(pooled,aes(Dose,frequency,colour=Profiler,linetype=Arm,shape=Arm,group=interaction(Profiler,Arm)))+
  geom_line(linewidth=MANUSCRIPT_LINE_MM)+geom_point(size=1.7)+scale_colour_profiler()+
  scale_linetype_manual(values=c(Observed='solid',Reference='dashed'))+
  scale_shape_manual(values=c(Observed=16,Reference=1))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
  labs(tag='A',x='Nominal per-target dose (%)',y='Positive significance frequency',
       subtitle='All targets/cohorts; n = 20/group')+theme_manuscript()+
  theme(legend.box='vertical')
target_panel <- function(tool,tag) {
  z <- low[low$profiler==tool,];z$Target <- factor(z$target_label,levels=rev(sort(unique(low$target_label))))
  segments <- reshape(z[,c('target_label','arm','frequency')],idvar='target_label',timevar='arm',direction='wide')
  segments$Target <- factor(segments$target_label,levels=levels(z$Target))
  ggplot(z,aes(frequency,Target,shape=Arm))+
    geom_segment(data=segments,aes(x=frequency.observed,xend=frequency.reference,y=Target,yend=Target),
                 inherit.aes=FALSE,colour='grey65',linewidth=MANUSCRIPT_AXIS_MM)+
    geom_point(colour=MANUSCRIPT_PROFILER_COLORS[tools[tool]],size=1.7)+
    scale_shape_manual(values=c(Observed=16,Reference=1))+
    scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+
    labs(tag=tag,x='Positive significance frequency',y=NULL,subtitle=tools[tool])+
    theme_manuscript()
}
b <- target_panel('kraken2_bracken','B');c <- target_panel('metaphlan4','C')
z <- ggplot(reason,aes(Profiler,frequency))+
  geom_col(fill='grey70',width=.55)+
  geom_text(aes(label=paste0(count,'/600\n',explanation)),vjust=-.25,size=7/ggplot2::.pt)+
  scale_y_continuous(limits=c(0,.25),breaks=c(0,.1,.2))+
  labs(tag='D',x=NULL,y='Non-estimable fraction',subtitle='Lowest dose: distinct NE reasons')+theme_manuscript()
panels <- list(a,b,c,z)
recolour <- function(g,fun) {
  if(!is.null(g$gp)) for(k in c('col','fill')) {
    v <- g$gp[[k]]
    if(is.character(v)) {ok <- !is.na(v)&v!='transparent';if(any(ok))v[ok]<-fun(v[ok]);g$gp[[k]]<-v}
  }
  if(!is.null(g$grobs))g$grobs<-lapply(g$grobs,recolour,fun=fun)
  if(!is.null(g$children))g$children<-do.call(grid::gList,lapply(g$children,recolour,fun=fun))
  g
}
draw <- function(transform=NULL) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2)))
  for(i in seq_along(panels)) {
    vp <- grid::viewport(layout.pos.row=(i-1)%/%2+1,layout.pos.col=(i-1)%%2+1)
    if(is.null(transform)) print(panels[[i]],vp=vp) else {
      grid::pushViewport(vp);grid::grid.draw(recolour(ggplotGrob(panels[[i]]),transform));grid::popViewport()
    }
  }
  grid::popViewport()
}
base <- file.path(out,'reference_response_compact')
cairo_pdf(paste0(base,'.pdf'),width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT);draw();dev.off()
svg(paste0(base,'.svg'),width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT);draw();dev.off()
png(paste0(base,'.png'),width=180,height=210,units='mm',res=450,type='cairo',bg='white');draw();dev.off()
stopifnot(requireNamespace('colorspace',quietly=TRUE))
cairo_pdf(file.path(out,'accessibility_review.pdf'),width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT)
for(fun in list(colorspace::desaturate,colorspace::deutan,colorspace::protan))draw(fun)
dev.off()
writeLines(c(
  'Pilot reference-response comparison in adenoma backgrounds, n=20 people per artificial group.',
  'A: positive coefficient and full-family MaAsLin2 q<=0.05 in observed versus baseline-anchored reference responses.',
  'Each point pools ten targets, three cohorts and 20 saved allocations: 600 target-context rows, not independent studies.',
  'B-C: lowest nominal per-target dose (0.001%); 60 conditional allocations/target, pooled across cohorts.',
  'D: observed NE states at that dose, each denominator 600; Bracken perfect fits are not absent signals.',
  'All significance denominators include NE; no confidence interval is assigned to repeated-draw frequencies.',
  'Reference is profiler-specific and baseline-biased, not common biological truth; no causal attribution to clinical non-association.',
  'Panel D intentionally uses a different frequency range to display small fractions; all bars begin at zero.',
  'Lines in A join tested doses; detailed cohort/species pages remain in the complete review booklet.',
  'Primary inference only; HC3 sensitivity and non-estimable states remain in full source reports.'),file.path(out,'CAPTION.txt'))
file.copy(style,file.path(out,'provenance','figure_style.R'))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
file.copy(script,file.path(out,'provenance','plot_reference_response_compact.R'))
writeLines(capture.output(sessionInfo()),file.path(out,'provenance','sessionInfo.txt'))
hash <- function(p) strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' ')[[1]][1]
writeLines(c(paste('reference_report_manifest',hash(file.path(report,'SHA256SUMS')),sep='\t'),
             paste('ne_audit_manifest',hash(file.path(audit,'SHA256SUMS')),sep='\t')),
           file.path(out,'provenance','input_hashes.tsv'))
files <- sort(list.files(out,recursive=TRUE,full.names=TRUE))
writeLines(vapply(files,function(p) paste0(hash(p),'  ',substring(p,nchar(out)+2)),character(1)),file.path(out,'SHA256SUMS'))
