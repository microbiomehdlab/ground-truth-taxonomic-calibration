#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop('Usage: task mode repository')
task <- args[1]; mode <- args[2]; repo <- args[3]
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/paired_difference_context.R'))
source(file.path(repo,'analysis_v2/lib/paired_null.R'))
source(file.path(repo,'analysis_v2/lib/unpaired_null.R'))
save_table <- function(x,name) write.table(x,file.path(task,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
if (mode=='sensitivity') {
  a <- as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE))
  storage.mode(a) <- 'double'
  meta <- read.delim(file.path(task,'metadata.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
  all <- lapply(c(1e-9,1e-8,1e-7),function(p) {
    x <- paired_difference_results(a,meta,colnames(a),p)$results
    x$pseudocount <- p; x
  })
  save_table(do.call(rbind,all),'sensitivity_results.tsv')
  # Preserve all features, not merely significant targets.
  summary <- do.call(rbind,lapply(seq_along(all),function(i)data.frame(
    pseudocount=c(1e-9,1e-8,1e-7)[i],features=nrow(all[[i]]),
    estimable=sum(all[[i]]$estimable),positive_discoveries=sum(all[[i]]$positive_discovery))))
  save_table(summary,'summary.tsv')
} else if (mode=='synthetic') {
  settings <- read.delim(file.path(task,'settings.tsv'))
  n <- settings$n_pairs; m <- settings$family_n; scenario <- settings$scenario
  set.seed(20261002L+n+match(scenario,c('gaussian','correlated_gaussian','sparse_symmetric','skewed_zero_mean'))*1000L)
  rows <- lapply(seq_len(1000),function(i) {
    d <- switch(scenario,
      gaussian=matrix(rnorm(n*m),n),
      correlated_gaussian=sqrt(.8)*matrix(rnorm(n*m),n)+sqrt(.2)*rnorm(n),
      sparse_symmetric=matrix(rnorm(n*m)*rbinom(n*m,1,.1),n),
      skewed_zero_mean=matrix(rexp(n*m)-1,n),stop('Unknown scenario'))
    stats <- paired_null_statistics(d,100*.Machine$double.eps*pmax(1,apply(abs(d),2,max)))
    data.frame(repetition=i,t(stats),check.names=FALSE)
  })
  draws <- do.call(rbind,rows); save_table(draws,'null_draws.tsv')
  count <- sum(draws$discoveries>0); ci <- binom.test(count,1000)$conf.int
  save_table(data.frame(scenario=scenario,n_pairs=n,family_n=m,repetitions=1000,
    any_discovery_rate=count/1000,mc_ci_low=ci[1],mc_ci_high=ci[2],
    mean_discoveries=mean(draws$discoveries),
    review_flag=count/1000>.075 || ci[1]>.05),'summary.tsv')
} else if (mode=='baseline_null') {
  source(file.path(repo,'analysis_v2/lib/maaslin_context.R'))
  if (!requireNamespace('Maaslin2',quietly=TRUE) || as.character(packageVersion('Maaslin2'))!='1.18.0')
    stop('Pinned MaAsLin2 1.18.0 required')
  a <- as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE)); storage.mode(a)<-'double'
  allocations <- read.delim(file.path(task,'allocations.tsv'),stringsAsFactors=FALSE)
  y <- log2(1+a/1e-8)
  keys <- unique(allocations[,c('pool','allocation_id','n')])
  draws <- vector('list',nrow(keys)); checked <- character()
  for (i in seq_len(nrow(keys))) {
    k <- keys[i,]; rows <- allocations[allocations$pool==k$pool & allocations$allocation_id==k$allocation_id & allocations$n==k$n,]
    if (anyDuplicated(rows$sample_id) || nrow(rows)!=2*k$n || !all(rows$sample_id %in% rownames(a))) stop('Invalid allocation')
    group <- as.integer(rows$group=='cases')
    if (sum(group)!=k$n || sum(group==0)!=k$n) stop('Invalid group counts')
    fit <- unpaired_null_results(y[rows$sample_id,,drop=FALSE],group)
    stratum <- paste(k$pool,k$n,sep='_')
    # Verify the fast LM calculation against the actual pinned primary wrapper
    # once for each n/pool stratum before trusting its repeated-null summaries.
    if (!stratum %in% checked) {
      meta <- data.frame(group=group,row.names=rows$sample_id)
      native <- fit_maaslin_context(a[rows$sample_id,,drop=FALSE],meta,colnames(a),file.path(task,paste0('backend_',stratum)))
      if (!identical(native$status,fit$status) || !isTRUE(all.equal(native$raw_p,fit$raw_p,tolerance=1e-7)) ||
          !isTRUE(all.equal(native$beta,fit$beta,tolerance=1e-7))) stop('Pinned backend/fast LM disagreement')
      checked <- c(checked,stratum)
    }
    draws[[i]] <- data.frame(pool=k$pool,allocation_id=k$allocation_id,n=k$n,
      discoveries=sum(fit$estimable & fit$wrapper_q<=.05),
      positive_discoveries=sum(fit$estimable & fit$wrapper_q<=.05 & !is.na(fit$beta) & fit$beta>0),
      estimable=sum(fit$estimable))
  }
  draws <- do.call(rbind,draws); save_table(draws,'null_draws.tsv')
  strata <- unique(draws[,c('pool','n')])
  summary <- do.call(rbind,lapply(seq_len(nrow(strata)),function(i) {
    s <- strata[i,]; x <- draws[draws$pool==s$pool & draws$n==s$n,]
    count <- sum(x$discoveries>0); rate <- count/nrow(x)
    ci <- if (s$pool=='full') binom.test(count,nrow(x))$conf.int else c(NA_real_,NA_real_)
    data.frame(pool=s$pool,n=s$n,repetitions=nrow(x),any_discovery_rate=rate,
      mc_ci_low=ci[1],mc_ci_high=ci[2],mean_discoveries=mean(x$discoveries),
      review_flag=rate>.075 || (!is.na(ci[1]) && ci[1]>.05))
  }))
  save_table(summary,'summary.tsv')
} else stop('Unknown validation mode')
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
cat('[PASS] Validation computation:',mode,'\n')
