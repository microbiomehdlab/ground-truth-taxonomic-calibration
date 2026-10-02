# Group-only equal-variance LM inference, checked against the pinned wrapper.
unpaired_null_results <- function(y,group) {
  if (!is.matrix(y) || any(!is.finite(y)) || length(group)!=nrow(y) ||
      !all(group %in% 0:1) || any(table(factor(group,levels=0:1))<2)) stop('Invalid unpaired design')
  control <- y[group==0,,drop=FALSE]; cases <- y[group==1,,drop=FALSE]
  beta <- colMeans(cases)-colMeans(control)
  residual <- rbind(sweep(control,2,colMeans(control),'-'),sweep(cases,2,colMeans(cases),'-'))
  tolerance <- 100*.Machine$double.eps*pmax(1,apply(abs(y),2,max))
  variable <- apply(y,2,function(v)diff(range(v)))>tolerance
  ok <- variable & apply(abs(residual),2,max)>tolerance
  se <- sqrt(colSums(residual^2)/(nrow(y)-2)*(1/nrow(control)+1/nrow(cases)))
  p <- rep(NA_real_,ncol(y)); p[ok] <- 2*pt(-abs(beta[ok]/se[ok]),nrow(y)-2)
  if (any(!is.finite(p[ok]))) stop('Invalid LM inference')
  status <- ifelse(!variable,'NON_ESTIMABLE_CONSTANT',ifelse(ok,'ESTIMABLE','NON_ESTIMABLE_PERFECT_FIXED_FIT'))
  bookkeeping <- ifelse(ok,p,1)
  data.frame(feature=colnames(y),beta=ifelse(ok,beta,NA_real_),raw_p=p,estimable=ok,
    status=status,wrapper_q=p.adjust(bookkeeping,'BH'),stringsAsFactors=FALSE)
}
