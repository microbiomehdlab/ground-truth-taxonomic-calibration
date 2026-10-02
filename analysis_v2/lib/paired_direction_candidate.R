# Different estimand from mean change: direction among nonzero changes.
# No symmetry requirement, but independent biological pairs and the sign-null
# Pr(D>0 | D!=0)=.5 are required. This does not certify a skewed mean-null test.
paired_direction_candidate <- function(d,tolerance) {
  if(!is.matrix(d) || nrow(d)<2 || any(!is.finite(d)) ||
     length(tolerance)!=ncol(d) || any(!is.finite(tolerance)) || any(tolerance<0)) stop('Invalid paired direction inputs')
  positive<-colSums(sweep(d,2,tolerance,'>'))
  negative<-colSums(sweep(d,2,-tolerance,'<'))
  n<-positive+negative
  raw<-rep(NA_real_,ncol(d));ok<-n>0
  raw[ok]<-pmin(1,2*pbinom(pmin(positive[ok],negative[ok]),n[ok],.5))
  bookkeeping<-ifelse(ok,raw,1)
  data.frame(feature=colnames(d),n_pairs=nrow(d),positive=positive,negative=negative,
    zero_or_numerical_tie=nrow(d)-n,informative_pairs=n,raw_p=raw,
    estimable=ok,full_family_bh_q=p.adjust(bookkeeping,'BH'),
    positive_discovery=ok & positive>negative & p.adjust(bookkeeping,'BH')<=.05,
    estimand='Pr(positive | nonzero change)',stringsAsFactors=FALSE)
}
