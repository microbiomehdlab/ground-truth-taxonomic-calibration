# Conditional sign-flip diagnostic, not empirical biological-null truth.
paired_null_statistics <- function(d, tolerance) {
  if (!is.matrix(d) || nrow(d)<2 || any(!is.finite(d)) ||
      length(tolerance)!=ncol(d) || any(!is.finite(tolerance)) || any(tolerance<0))
    stop('Invalid difference matrix or tolerances')
  n <- nrow(d); beta <- colMeans(d)
  se <- sqrt(colSums(sweep(d,2,beta,'-')^2)/(n-1)/n)
  ok <- se>tolerance & se>=10*.Machine$double.eps*abs(beta)
  p <- rep(1,ncol(d))
  p[ok] <- 2*pt(-abs(beta[ok]/se[ok]),df=n-1)
  q <- p.adjust(p,'BH')
  c(estimable=sum(ok),discoveries=sum(ok & q<=0.05),
    positive_discoveries=sum(ok & q<=0.05 & beta>0),
    nominal_rejections=sum(ok & p<=0.05))
}

paired_signflip_diagnostic <- function(d,tolerance,repetitions=1000L,seed=20261002L) {
  if (repetitions<1 || repetitions!=as.integer(repetitions)) stop('Invalid repetitions')
  set.seed(seed)
  # Same person-level sign vector across every species preserves dependence.
  signs <- matrix(sample(c(-1,1),nrow(d)*repetitions,replace=TRUE),nrow=nrow(d))
  draws <- t(vapply(seq_len(repetitions),function(i)
    paired_null_statistics(d*signs[,i],tolerance),numeric(4)))
  list(draws=data.frame(repetition=seq_len(repetitions),draws),signs=signs)
}
