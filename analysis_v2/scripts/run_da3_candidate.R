#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop('Usage: task mode repo')
task<-args[1];mode<-args[2];repo<-args[3]
source(file.path(repo,'analysis_v2/lib/exact_da3_candidate.R'))
source(file.path(repo,'analysis_v2/lib/unpaired_null.R'))
save_table<-function(x,name)write.table(x,file.path(task,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE));storage.mode(a)<-'double'
if(any(!is.finite(a)) || any(a<0 | a>1.00001))stop('Invalid native input')
y<-log2(1+a/1e-8)
evaluate<-function(selected,group){
  candidate<-exact_da3_candidate(selected,group)
  old<-unpaired_null_results(selected,group)
  list(features=candidate,counts=data.frame(
    parametric_bh=sum(old$estimable & old$wrapper_q<=.05),
    exact_bh=sum(candidate$variable & candidate$exact_bh_q<=.05),
    exact_max_fwer=sum(candidate$variable & candidate$max_statistic_fwer_p<=.05)))
}
if(mode=='positive'){
  meta<-read.delim(file.path(task,'metadata.tsv'),row.names=1,stringsAsFactors=FALSE)
  if(!setequal(rownames(meta),rownames(y)))stop('Metadata mismatch')
  meta<-meta[rownames(y),,drop=FALSE]
  if(anyDuplicated(meta$biological_sample_id))stop('DA3 requires distinct people')
  fit<-evaluate(y,as.integer(meta$group))
  save_table(fit$features,'candidate_results.tsv');save_table(fit$counts,'summary.tsv')
}else if(mode=='null'){
  allocations<-read.delim(file.path(task,'allocations.tsv'),stringsAsFactors=FALSE)
  keys<-unique(allocations[,c('pool','allocation_id')]);draws<-vector('list',nrow(keys))
  for(i in seq_len(nrow(keys))){
    key<-keys[i,];r<-allocations[allocations$pool==key$pool & allocations$allocation_id==key$allocation_id,]
    if(nrow(r)!=10 || anyDuplicated(r$sample_id) || !all(r$sample_id %in% rownames(y)))stop('Invalid allocation')
    fit<-evaluate(y[r$sample_id,,drop=FALSE],as.integer(r$group=='cases'))
    draws[[i]]<-cbind(key,fit$counts)
    if(i%%100==0)cat('[PROGRESS]',i,'/',nrow(keys),'allocations\n')
  }
  draws<-do.call(rbind,draws);save_table(draws,'null_draws.tsv')
  summary<-do.call(rbind,lapply(unique(draws$pool),function(pool){
    x<-draws[draws$pool==pool,];do.call(rbind,lapply(c('parametric_bh','exact_bh','exact_max_fwer'),function(method){
      count<-sum(x[[method]]>0);ci<-if(pool=='full')binom.test(count,nrow(x))$conf.int else c(NA_real_,NA_real_)
      data.frame(pool=pool,method=method,repetitions=nrow(x),any_discovery_rate=count/nrow(x),
        mean_discoveries=mean(x[[method]]),mc_ci_low=ci[1],mc_ci_high=ci[2])
    }))
  }))
  save_table(summary,'summary.tsv')
}else stop('Unknown mode')
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
cat('[PASS] Exact DA3 candidate computation:',mode,'\n')
