args <- commandArgs(TRUE); stopifnot(length(args)==1)
root <- args[1]
d <- read.delim(file.path(root,'dose_response_summary.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
palette <- c(yachida='#0072B2',feng='#D55E00',zeller='#009E73')
tools <- c(kraken2_bracken='Kraken2 / Bracken',metaphlan4='MetaPhlAn4')
metrics <- c(native_increment='Native fraction: spiked minus original',
             recovery_ratio='Reference-adjusted recovery ratio',legacy_reference_ratio='Legacy read-reference ratio')
for (pop in sort(unique(d$population))) for (metric in names(metrics)) {
  pdf(file.path(root,paste0(pop,'_',metric,'.pdf')),width=9,height=10,useDingbats=FALSE)
  for (label in sort(unique(d$target_label[d$population==pop]))) {
    par(mfrow=c(3,2),mar=c(4,4.6,2.6,1),oma=c(2,0,3,0),las=1,cex=.85)
    for (condition in c('Control','Adenoma','CRC')) for (tool in names(tools)) {
      z <- d[d$population==pop & d$target_label==label & d$condition==condition & d$profiler==tool,]
      stopifnot(nrow(z)>0)
      med <- paste0(metric,'_median'); low <- paste0(metric,'_q1'); high <- paste0(metric,'_q3')
      reference <- if(metric=='native_increment') 0 else 1
      ylim <- range(c(z[[low]],z[[high]],reference),finite=TRUE)
      if(diff(ylim)==0) ylim <- ylim+c(-.05,.05)
      plot(NA,xlim=range(z$nominal_target_dose),ylim=ylim,log='x',
           xlab='Nominal target read fraction',ylab=metrics[[metric]],
           main=paste(tools[[tool]],condition,sep=' | '),bty='l')
      abline(h=reference,col='grey65',lty=2)
      for (cohort in names(palette)) {
        v <- z[z$cohort==cohort,]; v <- v[order(v$nominal_target_dose),]
        lines(v$nominal_target_dose,v[[med]],col=palette[cohort],lwd=1.5,type='b',pch=16,cex=.65)
        segments(v$nominal_target_dose,v[[low]],v$nominal_target_dose,v[[high]],col=palette[cohort])
      }
      if(condition=='Control' && tool==names(tools)[1]) legend('topleft',legend=names(palette),
          col=palette,lty=1,pch=16,bty='n',cex=.75)
    }
    mtext(paste(label,pop,'paired dose response',sep=' | '),outer=TRUE,line=1,font=2)
    mtext('Median and person IQR; not confidence intervals. Reference scales differ by profiler.',outer=TRUE,side=1,cex=.8)
  }
  dev.off()
}
cat('[PASS] Six DA2 multipage PDFs; complete plotted values retained in dose_response_summary.tsv\n')
