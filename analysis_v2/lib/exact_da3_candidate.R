# Candidate for 5+5 randomized-group experiments only, not clinical DA1.
# Two-sided absolute mean contrasts, scaled by label-invariant pooled SD.
exact_da3_candidate <- function(y,group) {
  if (!is.matrix(y) || nrow(y)!=10 || any(!is.finite(y)) ||
      length(group)!=10 || anyNA(group) || !all(group %in% 0:1) || sum(group)!=5 ||
      is.null(colnames(y)) || anyDuplicated(colnames(y))) stop('Exact candidate requires complete 5+5 design')
  combinations <- combn(10,5)
  weights <- matrix(-.2,nrow=252,ncol=10)
  for (i in seq_len(252)) weights[i,combinations[,i]] <- .2
  centered <- sweep(y,2,colMeans(y),'-')
  scale <- sqrt(colSums(centered^2)/9)
  tolerance <- 100*.Machine$double.eps*pmax(1,apply(abs(y),2,max))
  variable <- scale>tolerance
  standardized <- centered
  standardized[,variable] <- sweep(centered[,variable,drop=FALSE],2,scale[variable],'/')
  standardized[,!variable] <- 0
  null <- abs(weights %*% standardized)
  observed <- abs(drop((ifelse(group==1,.2,-.2)) %*% standardized))
  epsilon <- 100*.Machine$double.eps*pmax(1,observed)
  p <- colMeans(sweep(null,2,observed-epsilon,'>=') )
  maximum <- apply(null,1,max)
  adjusted <- vapply(seq_along(observed),function(j)mean(maximum>=observed[j]-epsilon[j]),numeric(1))
  p[!variable] <- 1; adjusted[!variable] <- 1
  beta <- colMeans(y[group==1,,drop=FALSE])-colMeans(y[group==0,,drop=FALSE])
  data.frame(feature=colnames(y),beta=beta,variable=variable,exact_p=p,
    exact_bh_q=p.adjust(p,'BH'),max_statistic_fwer_p=adjusted,
    permutations=252,stringsAsFactors=FALSE)
}
