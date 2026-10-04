#!/usr/bin/env Rscript
fit_clinical_spike <- function(root,repo,backend=NULL) {
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/maaslin_context.R'))
a <- read.delim(file.path(root,'abundance.tsv'),check.names=FALSE,row.names=1)
m <- read.delim(file.path(root,'metadata.tsv'),check.names=FALSE,row.names=1)
m$sex <- factor(m$sex,levels=c('Female','Male'))
stopifnot(!anyNA(m$sex), all(m$spike_state == ifelse(m$group==1,'spiked','original')))
a <- as.matrix(a)
r <- fit_maaslin_context(a,m,colnames(a),file.path(root,'fit'),
                        paired=FALSE,covariates=c('age','sex'),backend=backend)
m <- m[rownames(a),,drop=FALSE]
design <- model.matrix(~group+age+sex,m)
df <- nrow(design)-ncol(design)
# Independent arithmetic audit of every estimable native package result.
transformed <- log2(1+a/1e-8)
for (i in which(r$estimable)) {
  fit <- lm.fit(design,transformed[,i])
  inverse <- chol2inv(qr.R(fit$qr))
  se <- sqrt(sum(fit$residuals^2)/df * inverse[2,2])
  beta <- unname(fit$coefficients[2])
  p <- 2*pt(abs(beta/se),df,lower.tail=FALSE)
  observed <- c(r$beta[i],r$stderr[i],r$raw_p[i])
  expected <- c(beta,se,p)
  if (any(abs(observed-expected)>1e-10+1e-7*pmax(abs(observed),abs(expected))))
    stop('Independent LM audit failed: ',r$feature[i])
}
r$residual_df <- df
r$ci95_low <- r$beta-qt(.975,df)*r$stderr
r$ci95_high <- r$beta+qt(.975,df)*r$stderr
r$significant_two_sided <- r$estimable & r$wrapper_q<=.05
for (group in 0:1) {
  x <- a[m$group==group,,drop=FALSE]
  prefix <- if(group==0) 'control' else 'adenoma_spiked'
  r[[paste0(prefix,'_n')]] <- nrow(x)
  r[[paste0(prefix,'_prevalence')]] <- colMeans(x>0)
  r[[paste0(prefix,'_native_mean')]] <- colMeans(x)
  r[[paste0(prefix,'_native_median')]] <- apply(x,2,median)
}
write.table(r,file.path(root,'clinical_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
invisible(r)
}
if (sys.nframe()==0) {
  args <- commandArgs(TRUE)
  stopifnot(length(args)==2)
  fit_clinical_spike(args[1],args[2])
}
