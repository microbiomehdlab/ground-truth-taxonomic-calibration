# Descriptive conditional allocations, no new inference or interval estimation.
args <- commandArgs(TRUE); stopifnot(length(args)==2)
root <- normalizePath(args[1]); source(normalizePath(args[2])); library(ggplot2)
f <- read.delim(file.path(root,'frequencies.tsv'),stringsAsFactors=FALSE)
e <- read.delim(file.path(root,'effects.tsv'),stringsAsFactors=FALSE)
s <- read.delim(file.path(root,'significance.tsv'),stringsAsFactors=FALSE)
tools <- c(kraken2_bracken='Kraken2/Bracken',metaphlan4='MetaPhlAn 4')
dir.create(file.path(root,'figures')); dir.create(file.path(root,'accessibility'))
if(!requireNamespace('colorspace',quietly=TRUE)) stop('colorspace required for accessibility review')
recolour <- function(g,fun) {
  if(!is.null(g$gp)) for(k in c('col','fill')) {
    v <- g$gp[[k]]
    if(is.character(v)) {ok <- !is.na(v) & v!='transparent';if(any(ok)) v[ok] <- fun(v[ok]);g$gp[[k]] <- v}
  }
  if(!is.null(g$grobs)) g$grobs <- lapply(g$grobs,recolour,fun=fun)
  if(!is.null(g$children)) g$children <- do.call(grid::gList,lapply(g$children,recolour,fun=fun))
  g
}
book <- file.path(root,'reference_response_review.pdf')
grDevices::cairo_pdf(book,width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT)
book_device <- dev.cur()
save <- function(p,name) {
  p <- p+theme_manuscript(); dev.set(book_device); print(p)
  base <- file.path(root,'figures',name)
  cairo_pdf(paste0(base,'.pdf'),width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT);print(p);dev.off()
  svg(paste0(base,'.svg'),width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT);print(p);dev.off()
  png(paste0(base,'.png'),width=180,height=210,units='mm',res=450,type='cairo',bg='white');print(p);dev.off()
  cairo_pdf(file.path(root,'accessibility',paste0(name,'.pdf')),width=180/25.4,height=210/25.4,family=MANUSCRIPT_FONT)
  for(fun in list(colorspace::deutan,colorspace::protan,colorspace::desaturate)) {
    grid::grid.newpage();grid::grid.draw(recolour(ggplotGrob(p),fun))
  };dev.off()
}
decorate <- function(d) {
  d$Dose <- factor(d$nominal_total_dose,levels=c(.0001,.001,.01),labels=c('0.001','0.01','0.1'))
  d$Target <- factor(d$target_label,levels=sort(unique(f$target_label)))
  d$N <- paste0('n = ',d$n,' / group');d
}
for(cohort in c('yachida','feng','zeller')) for(tool in names(tools)) {
  d <- decorate(f[f$cohort==cohort & f$profiler==tool & f$inference=='primary_MaAsLin2',])
  # Seven exhaustive categories remain in source tables. Significance and NE
  # are separate views, each denominator includes all 20 allocations.
  curves <- decorate(s[s$cohort==cohort & s$profiler==tool,])
  curves$Series <- ifelse(curves$arm=='observed','Observed','Reference')
  colour <- MANUSCRIPT_PROFILER_COLORS[tools[tool]]
  p <- ggplot(curves,aes(Dose,frequency,group=Series,linetype=Series,shape=Series))+
    geom_line(colour=colour,linewidth=MANUSCRIPT_LINE_MM)+geom_point(colour=colour,size=1.7)+
    facet_wrap(vars(Target,N),ncol=4)+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
    scale_linetype_manual(values=c(Observed='solid',Reference='dashed'))+
    scale_shape_manual(values=c(Observed=16,Reference=1))+
    labs(x='Nominal per-target read-pair dose (%)',y='Positive significance / all allocations',
      title=paste(cohort,tools[tool]),subtitle='MaAsLin2; full-family q ≤ 0.05; 20 conditional allocations/cell')
  save(p,paste(cohort,tool,'significance',sep='_'))
  ne <- d[d$category %in% c('OBSERVED_NON_ESTIMABLE','REFERENCE_NON_ESTIMABLE','BOTH_NON_ESTIMABLE'),]
  ne$State <- factor(ne$category,levels=c('OBSERVED_NON_ESTIMABLE','REFERENCE_NON_ESTIMABLE','BOTH_NON_ESTIMABLE'),
                    labels=c('Observed only NE','Reference only NE','Both NE'))
  p <- ggplot(ne,aes(Dose,frequency,shape=State,linetype=State,group=State))+
    geom_line(linewidth=MANUSCRIPT_LINE_MM)+geom_point(size=1.7)+facet_wrap(vars(Target,N),ncol=4)+
    scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
    labs(x='Nominal per-target read-pair dose (%)',y='Non-estimable / all allocations',
      title=paste(cohort,tools[tool]),subtitle='NE is not ordinary nonsignificance; 20 allocations/cell')
  save(p,paste(cohort,tool,'nonestimable',sep='_'))
  for(quantity in c('beta_difference','observed_stderr','reference_stderr')) {
    z <- decorate(e[e$cohort==cohort & e$profiler==tool & e$quantity==quantity,])
    p <- ggplot(z,aes(Dose,median,group=1))+geom_hline(yintercept=0,linetype='dashed',linewidth=MANUSCRIPT_AXIS_MM)+
      geom_linerange(aes(ymin=lower,ymax=upper),colour=colour,linewidth=MANUSCRIPT_AXIS_MM)+
      geom_line(colour=colour,linewidth=MANUSCRIPT_LINE_MM)+geom_point(colour=colour,size=1.7)+
      geom_text(aes(label=paste0('m=',matched_estimable)),y=Inf,vjust=1.3,size=7/ggplot2::.pt)+
      facet_wrap(vars(Target,N),ncol=4)+labs(x='Nominal per-target read-pair dose (%)',
        y=if(quantity=='beta_difference') 'Observed − reference coefficient' else paste(quantity,'(transformed units)'),
        title=paste(cohort,tools[tool]),subtitle='Median/IQR across matched estimable allocations; not confidence intervals')
    save(p,paste(cohort,tool,quantity,sep='_'))
  }
}
dev.set(book_device);dev.off()
writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
dir.create(file.path(root,'provenance'),showWarnings=FALSE)
file.copy(args[2],file.path(root,'provenance','figure_style.R'),overwrite=TRUE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
file.copy(script,file.path(root,'provenance','plot_reference_response.R'),overwrite=TRUE)
files <- sort(list.files(root,recursive=TRUE,full.names=TRUE))
files <- files[basename(files)!='SHA256SUMS']
lines <- vapply(files,function(p) paste0(strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' ')[[1]][1],
  '  ',substring(p,nchar(root)+2)),character(1))
writeLines(lines,file.path(root,'SHA256SUMS'))
