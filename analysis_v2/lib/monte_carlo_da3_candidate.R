# Balanced independent groups; exchangeable-label global-null candidate.
# Fixed resample budget and plus-one p values. Never call these exact p values.
monte_carlo_da3_candidate <- function(y,group,repetitions=9999L,seed=1L,chunk=256L) {
  if (!is.matrix(y) || any(!is.finite(y)) || length(group)!=nrow(y) || anyNA(group) ||
      !all(group %in% 0:1) || sum(group)!=sum(group==0) || sum(group)<2 ||
      is.null(colnames(y)) || anyDuplicated(colnames(y)) || repetitions<1 ||
      repetitions!=as.integer(repetitions) || chunk<1 || chunk!=as.integer(chunk)) stop('Invalid balanced permutation context')
  total<-nrow(y);n<-sum(group);centered<-sweep(y,2,colMeans(y),'-')
  scale<-sqrt(colSums(centered^2)/(total-1))
  tolerance<-100*.Machine$double.eps*pmax(1,apply(abs(y),2,max));variable<-scale>tolerance
  standardized<-centered;standardized[,variable]<-sweep(centered[,variable,drop=FALSE],2,scale[variable],'/')
  standardized[,!variable]<-0
  observed<-abs(drop(ifelse(group==1,1/n,-1/n)%*%standardized))
  epsilon<-100*.Machine$double.eps*pmax(1,observed)
  raw_count<-numeric(ncol(y));maximum_count<-numeric(ncol(y))
  set.seed(seed)
  for(start in seq.int(1L,repetitions,by=chunk)) {
    count<-min(chunk,repetitions-start+1L)
    weights<-matrix(-1/n,count,total)
    for(i in seq_len(count)) weights[i,sample.int(total,n)]<-1/n
    null<-abs(weights%*%standardized)
    raw_count<-raw_count+colSums(sweep(null,2,observed-epsilon,'>='))
    maximum<-apply(null,1,max)
    maximum_count<-maximum_count+vapply(seq_along(observed),function(j)sum(maximum>=observed[j]-epsilon[j]),numeric(1))
  }
  p<-(raw_count+1)/(repetitions+1);adjusted<-(maximum_count+1)/(repetitions+1)
  p[!variable]<-1;adjusted[!variable]<-1
  data.frame(feature=colnames(y),beta=colMeans(y[group==1,,drop=FALSE])-colMeans(y[group==0,,drop=FALSE]),
    variable=variable,mc_p=p,mc_bh_q=p.adjust(p,'BH'),max_statistic_fwer_p=adjusted,
    raw_exceedances=raw_count,max_exceedances=maximum_count,resamples=repetitions,
    minimum_mc_p=1/(repetitions+1),seed=seed,stringsAsFactors=FALSE)
}
