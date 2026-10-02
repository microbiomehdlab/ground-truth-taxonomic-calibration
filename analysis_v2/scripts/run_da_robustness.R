#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3)stop('Usage: task mode repo')
task<-args[1];mode<-args[2];repo<-args[3]
source(file.path(repo,'analysis_v2/lib/monte_carlo_da3_candidate.R'))
source(file.path(repo,'analysis_v2/lib/paired_direction_candidate.R'))
source(file.path(repo,'analysis_v2/lib/unpaired_null.R'))
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/paired_difference_context.R'))
save_table<-function(x,name)write.table(x,file.path(task,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
if(mode=='larger_null'){
  a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE));storage.mode(a)<-'double'
  if(any(!is.finite(a)) || any(a<0 | a>1.00001))stop('Invalid native input')
  y<-log2(1+a/1e-8);alloc<-read.delim(file.path(task,'allocations.tsv'),stringsAsFactors=FALSE)
  keys<-unique(alloc[,c('allocation_id','seed')]);draws<-vector('list',nrow(keys))
  for(i in seq_len(nrow(keys))){
    r<-alloc[alloc$allocation_id==keys$allocation_id[i],]
    if(anyDuplicated(r$sample_id) || !all(r$sample_id %in% rownames(y)) || !all(r$group %in% c('cases','controls')))stop('Invalid allocation')
    group<-as.integer(r$group=='cases');selected<-y[r$sample_id,,drop=FALSE]
    fit<-monte_carlo_da3_candidate(selected,group,seed=keys$seed[i])
    old<-unpaired_null_results(selected,group)
    draws[[i]]<-data.frame(allocation_id=keys$allocation_id[i],seed=keys$seed[i],n=sum(group),
      parametric_bh=sum(old$estimable & old$wrapper_q<=.05),
      mc_bh=sum(fit$variable & fit$mc_bh_q<=.05),
      mc_max_fwer=sum(fit$variable & fit$max_statistic_fwer_p<=.05),
      permutations=9999,minimum_mc_p=.0001)
    if(i==1)save_table(fit,'first_allocation_results.tsv')
    if(i%%10==0)cat('[PROGRESS]',i,'/',nrow(keys),'allocations\n')
  }
  save_table(do.call(rbind,draws),'draws.tsv')
}else if(mode=='paired_direction'){
  a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE));storage.mode(a)<-'double'
  meta<-read.delim(file.path(task,'metadata.tsv'),row.names=1,stringsAsFactors=FALSE)
  paired<-paired_difference_results(a,meta,colnames(a))
  direction<-paired_direction_candidate(paired$differences,paired$results$numerical_tolerance)
  save_table(direction,'direction_results.tsv')
  save_table(paired$results,'mean_results.tsv')
  save_table(data.frame(mean_positive_discoveries=sum(paired$results$positive_discovery),
    direction_positive_discoveries=sum(direction$positive_discovery)), 'summary.tsv')
}else if(mode=='direction_synthetic'){
  settings<-read.delim(file.path(task,'settings.tsv'));n<-settings$n;scenario<-settings$scenario;m<-3471
  set.seed(20261002L+n+match(scenario,c('gaussian','sparse_symmetric','skewed_sign_null','skewed_mean_null'))*1000L)
  draws<-lapply(seq_len(1000),function(i){
    d<-switch(scenario,gaussian=matrix(rnorm(n*m),n),
      sparse_symmetric=matrix(rnorm(n*m)*rbinom(n*m,1,.1),n),
      skewed_sign_null=matrix(rexp(n*m)-log(2),n),
      skewed_mean_null=matrix(rexp(n*m)-1,n),stop('Unknown scenario'))
    colnames(d)<-paste0('f',seq_len(m))
    fit<-paired_direction_candidate(d,100*.Machine$double.eps*pmax(1,apply(abs(d),2,max)))
    data.frame(repetition=i,discoveries=sum(fit$estimable & fit$full_family_bh_q<=.05))
  });draws<-do.call(rbind,draws);save_table(draws,'draws.tsv')
  count<-sum(draws$discoveries>0);ci<-binom.test(count,1000)$conf.int
  save_table(data.frame(n=n,scenario=scenario,repetitions=1000,any_discovery_rate=count/1000,
    mc_ci_low=ci[1],mc_ci_high=ci[2],
    interpretation=if(scenario=='skewed_mean_null')'NOT_A_SIGN_NULL_DIFFERENT_ESTIMAND' else 'SIGN_NULL_STRESS_TEST'),'summary.tsv')
}else stop('Unknown mode')
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
cat('[PASS] Robustness candidate:',mode,'\n')
