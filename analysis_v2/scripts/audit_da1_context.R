#!/usr/bin/env Rscript
# Audit an existing MaAsLin2 clinical fit; never replace its inference.
audit_da1_context <- function(task, out) {
  read <- function(name) read.delim(file.path(task,name),check.names=FALSE,
      stringsAsFactors=FALSE,quote='',comment.char='')
  close <- function(a,b,label) {
    if (length(a)!=length(b) || any(is.na(a)!=is.na(b)) ||
        any(abs(a-b)>1e-10+1e-7*pmax(abs(a),abs(b)),na.rm=TRUE))
      stop('Saved inference mismatch: ',label)
  }
  c <- read('context.tsv'); m <- read('metadata.tsv'); a <- read('abundance.tsv')
  r <- read('fit/context_results.tsv'); map <- read('fit/feature_map.tsv')
  stopifnot(nrow(c)==1,c$analysis=='DA1',c$paired==0,c$covariates=='age,sex',
      c$population=='community',as.numeric(c$nominal_total_dose)==0,
      !anyDuplicated(m$observation_id),!anyDuplicated(m$biological_sample_id),
      identical(a[[1]],m$observation_id),all(m$spike_state=='original'),
      setequal(m$group,0:1),all(m$sex %in% c('Female','Male')),
      all(is.finite(m$age)),all(m$age>0 & m$age<=120),
      identical(names(a)[-1],r$feature),identical(map$feature,r$feature),
      identical(map$model_feature,r$model_feature),!anyDuplicated(r$feature),
      !anyDuplicated(r$model_feature),nrow(m)==c$observations_n,nrow(r)==c$family_n)
  x <- as.matrix(a[,-1,drop=FALSE]); storage.mode(x) <- 'double'
  stopifnot(all(is.finite(x)),all(x>=0 & x<=1))
  m$sex <- factor(m$sex,levels=c('Female','Male'))
  design <- model.matrix(~group+age+sex,m)
  stopifnot(qr(design)$rank==4,nrow(m)>4)
  df <- nrow(m)-4
  stopifnot(!anyNA(r$estimable),!anyNA(r$status),
      all(is.finite(r$beta[r$estimable])),all(is.finite(r$stderr[r$estimable])),
      all(r$stderr[r$estimable]>0),all(is.finite(r$raw_p[r$estimable])),
      all(r$raw_p[r$estimable]>=0 & r$raw_p[r$estimable]<=1),
      all(is.finite(r$wrapper_q)),all(is.finite(r$p_for_BH)))
  native_path <- file.path(task,'fit/maaslin_native/all_results.tsv')
  native <- if(file.exists(native_path)) read('fit/maaslin_native/all_results.tsv') else NULL
  if (any(r$estimable)) {
    stopifnot(!is.null(native))
    native <- native[native$metadata=='group',,drop=FALSE]
    stopifnot(!anyDuplicated(native$feature),setequal(native$feature,r$model_feature[r$estimable]))
    j <- match(r$model_feature[r$estimable],native$feature)
    for (pair in list(c('beta','coef'),c('stderr','stderr'),c('raw_p','pval'),c('native_q','qval')))
      close(r[[pair[1]]][r$estimable],native[[pair[2]]][j],pair[1])
  }
  y <- log2(1+x/1e-8)
  for (i in seq_len(ncol(y))) {
    v <- y[,i]; tolerance <- 100*.Machine$double.eps*max(1,max(abs(v)))
    fit <- lm.fit(design,v)
    status <- if(diff(range(v))<=tolerance) 'NON_ESTIMABLE_CONSTANT' else
      if(max(abs(fit$residuals))<=tolerance) 'NON_ESTIMABLE_PERFECT_FIXED_FIT' else 'ESTIMABLE'
    stopifnot(r$status[i]==status,r$estimable[i]==(status=='ESTIMABLE'))
    if (status=='ESTIMABLE') {
      # Independently verify the frozen age/sex-adjusted OLS contract. Retain
      # saved MaAsLin2 numbers, not these recomputed numbers, in final tables.
      se <- sqrt(sum(fit$residuals^2)/df * solve(crossprod(design))[2,2])
      beta <- unname(fit$coefficients[2]); p <- 2*pt(-abs(beta/se),df)
      close(c(r$beta[i],r$stderr[i],r$raw_p[i]),c(beta,se,p),r$feature[i])
    } else stopifnot(all(is.na(unlist(r[i,c('beta','stderr','raw_p','native_q')]))))
  }
  bh <- ifelse(r$estimable,r$raw_p,1)
  close(r$p_for_BH,bh,'bookkeeping p'); close(r$wrapper_q,p.adjust(bh,'BH'),'full-family BH')
  stopifnot(identical(r$positive_discovery,r$estimable & !is.na(r$beta) & r$beta>0 & r$wrapper_q<=.05))
  r$residual_df <- df
  r$ci95_low <- r$beta-qt(.975,df)*r$stderr
  r$ci95_high <- r$beta+qt(.975,df)*r$stderr
  r$significant_two_sided <- r$estimable & r$wrapper_q<=.05
  r$effect_direction <- ifelse(!r$estimable,'not_estimable',ifelse(r$beta>0,'higher_in_disease',
      ifelse(r$beta<0,'lower_in_disease','zero')))
  for (g in 0:1) {
    prefix <- if(g==0) 'control' else 'disease'
    z <- x[m$group==g,,drop=FALSE]
    r[[paste0(prefix,'_n')]] <- nrow(z)
    r[[paste0(prefix,'_positive_n')]] <- colSums(z>0)
    r[[paste0(prefix,'_prevalence')]] <- colMeans(z>0)
    r[[paste0(prefix,'_native_mean')]] <- colMeans(z)
    r[[paste0(prefix,'_native_median')]] <- apply(z,2,median)
  }
  write.table(r,out,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  cat('[PASS] Saved DA1 inference and age/sex-adjusted design:',c$context_id,'\n')
}
if (sys.nframe()==0) {
  args <- commandArgs(TRUE); stopifnot(length(args)==2)
  audit_da1_context(args[1],args[2])
}
