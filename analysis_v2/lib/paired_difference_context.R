# DA2 primary: biological-person differences; no random-intercept optimizer.
paired_difference_results <- function(abundance, metadata, family, pseudocount=1e-8) {
  input <- prepare_maaslin_input(abundance,metadata,family,paired=TRUE,pseudocount=pseudocount)
  meta <- input$metadata
  if (!'group' %in% names(meta) || anyNA(meta$group) ||
      !all(as.character(meta$group)==ifelse(meta$spike_state=='original','0','1')))
    stop('Group must be 0 original and 1 spiked')
  people <- sort(unique(as.character(meta$biological_sample_id)))
  if (length(people)<2) stop('At least two biological pairs required')
  keys <- paste(meta$biological_sample_id,meta$spike_state,sep='\r')
  original <- match(paste(people,'original',sep='\r'),keys)
  spiked <- match(paste(people,'spiked',sep='\r'),keys)
  y <- as.matrix(input$data)
  differences <- y[spiked,,drop=FALSE]-y[original,,drop=FALSE]
  rownames(differences) <- people
  rows <- lapply(seq_along(family),function(j) {
    d <- differences[,j]; n <- length(d); beta <- mean(d)
    se <- sd(d)/sqrt(n)
    # Conservative numerical guard in transformed-response units, including
    # subtraction precision; never assign a fabricated p for constant changes.
    tolerance <- 100*.Machine$double.eps*max(1,abs(y[,j]))
    ok <- is.finite(se) && se>tolerance && se>=10*.Machine$double.eps*abs(beta)
    if (ok) {
      fit <- t.test(d,mu=0,alternative='two.sided',conf.level=0.95)
      if (!is.finite(fit$p.value)) stop('Invalid paired p value')
      p <- fit$p.value; ci <- fit$conf.int
      status <- 'ESTIMABLE'
    } else {
      p <- NA_real_; ci <- c(NA_real_,NA_real_)
      status <- if (max(abs(d))<=tolerance) 'NON_ESTIMABLE_NO_CHANGE' else 'NON_ESTIMABLE_CONSTANT_DIFFERENCE'
    }
    data.frame(feature=family[j],beta=beta,stderr=se,n_pairs=n,df=n-1,
      difference_sd=sd(d),numerical_tolerance=tolerance,ci_low=ci[1],ci_high=ci[2],
      raw_p=p,native_q=NA_real_,estimable=ok,status=status,stringsAsFactors=FALSE)
  })
  result <- do.call(rbind,rows)
  bh <- full_family_bh(family,result$feature,result$raw_p,result$estimable)
  result$p_for_BH <- bh$p_for_BH
  result$wrapper_q <- bh$wrapper_q
  result$positive_discovery <- result$estimable & result$beta>0 & result$wrapper_q<=0.05
  list(results=result,differences=differences)
}
