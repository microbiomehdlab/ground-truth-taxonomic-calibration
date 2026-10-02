#!/usr/bin/env Rscript
source('analysis_v2/lib/maaslin_contract.R')
source('analysis_v2/lib/paired_difference_context.R')
ids <- paste0('o',1:12)
meta <- data.frame(biological_sample_id=rep(paste0('person',1:6),2),
  spike_state=rep(c('original','spiked'),each=6),group=rep(0:1,each=6),row.names=ids)
y <- cbind(positive=c(1:6,1:6+c(1,2,1,3,2,4)),negative=c(7:12,7:12-c(1,2,1,3,2,4)),
  unchanged=rep(2,12),constant=c(rep(2,6),rep(3,6)))
a <- (2^y-1)*1e-8; rownames(a) <- ids
family <- colnames(a)
x <- paired_difference_results(a,meta,family)
stopifnot(nrow(x$results)==4,all(x$results$n_pairs==6),all(x$results$df==5))
for (j in 1:2) {
  reference <- t.test(y[7:12,j]-y[1:6,j])
  stopifnot(isTRUE(all.equal(x$results$raw_p[j],reference$p.value)),
    isTRUE(all.equal(x$results$ci_low[j],unname(reference$conf.int[1]))))
}
stopifnot(all(is.na(x$results$raw_p[3:4])),all(x$results$p_for_BH[3:4]==1),
  identical(x$results$status[3:4],c('NON_ESTIMABLE_NO_CHANGE','NON_ESTIMABLE_CONSTANT_DIFFERENCE')),
  isTRUE(all.equal(x$results$wrapper_q,p.adjust(x$results$p_for_BH,'BH'))))
order <- c(12:7,6:1)
reordered <- paired_difference_results(a[order,],meta[rev(order),],rev(family))
stopifnot(isTRUE(all.equal(x$results$raw_p,rev(reordered$results$raw_p))))
must_fail <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=function(e)e),'error'))
bad <- meta; bad$group[1] <- 1
must_fail(paired_difference_results(a,bad,family))
bad <- meta; bad$biological_sample_id[1] <- 'other'
must_fail(paired_difference_results(a,bad,family))
bad <- a; bad[1,1] <- NA_real_
must_fail(paired_difference_results(bad,meta,family))
must_fail(paired_difference_results(a[1:6,],meta[1:6,],family))
task <- tempfile('paired-runner-'); dir.create(task)
write.table(data.frame(analysis='DA2',paired=1),file.path(task,'context.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(observation_id=ids,a,check.names=FALSE),file.path(task,'abundance.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(observation_id=ids,meta),file.path(task,'metadata.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
code <- system2(file.path(R.home('bin'),'Rscript'),c('analysis_v2/scripts/run_paired_difference_context.R',shQuote(task),shQuote(getwd())))
stopifnot(code==0,file.exists(file.path(task,'fit/SUCCESS')))
saved <- read.delim(file.path(task,'fit/context_results.tsv'),check.names=FALSE)
stopifnot(nrow(saved)==4,isTRUE(all.equal(saved$raw_p,x$results$raw_p)))
stopifnot(file.exists(file.path(task,'fit/paired_differences.tsv')),file.exists(file.path(task,'sessionInfo.txt')))
unlink(task,recursive=TRUE)
cat('[PASS] Paired matching, t-test equivalence, full-family BH, degeneracy and invalid-input checks\n')
