# Source maaslin_contract.R before this file. Development wrapper, not frozen.
fit_maaslin_context <- function(abundance, metadata, family, output,
                               paired=FALSE, covariates=character(),
                               pseudocount=1e-8, backend=NULL) {
  if (file.exists(output)) stop('Choose a fresh output directory')
  input <- prepare_maaslin_input(abundance,metadata,family,paired,pseudocount)
  meta <- input$metadata
  if (!'group' %in% names(meta) || anyNA(meta$group) ||
      !all(meta$group %in% c(0,1)) || length(unique(meta$group))!=2)
    stop('group must encode reference=0 and intervention=1')
  if (paired && any(meta$group != as.integer(meta$spike_state=='spiked')))
    stop('Paired group disagrees with spike_state')
  if (anyDuplicated(covariates) || any(covariates %in% c('group','biological_sample_id')) ||
      !all(covariates %in% names(meta))) stop('Invalid covariates')
  variables <- c('group',covariates)
  if (anyNA(meta[,variables,drop=FALSE])) stop('Missing model metadata')
  design <- model.matrix(reformulate(variables),meta)
  if (qr(design)$rank!=ncol(design) || nrow(design)<=ncol(design))
    stop('Rank-deficient or saturated fixed-effect design')
  if (is.null(backend)) {
    if (!requireNamespace('Maaslin2',quietly=TRUE)) stop('MaAsLin2 not installed')
    if (as.character(packageVersion('Maaslin2'))!='1.18.0') stop('Require MaAsLin2 1.18.0')
    backend <- Maaslin2::Maaslin2
  }
  # Stable neutral IDs avoid punctuation/name repairs merging biological taxa.
  ids <- sprintf('feature_%06d',seq_along(family))
  data <- input$data
  names(data) <- ids
  results <- data.frame(feature=family,model_feature=ids,beta=NA_real_,stderr=NA_real_,
                        raw_p=NA_real_,native_q=NA_real_,estimable=FALSE,
                        status='NON_ESTIMABLE_CONSTANT',stringsAsFactors=FALSE)
  variable <- vapply(data,function(y) diff(range(y)) >
                      100*.Machine$double.eps*max(1,max(abs(y))),logical(1))
  # Pre-screen deterministic perfect fixed fits; never create machine-min p.
  for (i in which(variable)) {
    y <- data[[i]]
    residual <- lm.fit(design,y)$residuals
    if (max(abs(residual)) <= 100*.Machine$double.eps*max(1,max(abs(y)))) {
      variable[i] <- FALSE
      results$status[i] <- 'NON_ESTIMABLE_PERFECT_FIXED_FIT'
    }
  }
  dir.create(output,recursive=TRUE)
  write.table(data.frame(model_feature=ids,feature=family),file.path(output,'feature_map.tsv'),
              sep='\t',row.names=FALSE,quote=FALSE)
  if (any(variable)) {
    warnings <- character()
    fitout <- file.path(output,'maaslin_native')
    withCallingHandlers(backend(
      input_data=data[,variable,drop=FALSE],input_metadata=meta,output=fitout,
      fixed_effects=variables,
      random_effects=if(paired) 'biological_sample_id' else NULL,
      normalization='NONE',transform='NONE',analysis_method='LM',
      min_abundance=0,min_prevalence=0,min_variance=0,standardize=FALSE,
      correction='BH',max_significance=.05,cores=1,
      plot_heatmap=FALSE,plot_scatter=FALSE,save_models=TRUE),
      warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning') })
    if (length(warnings)) {
      writeLines(warnings,file.path(output,'warnings.txt'))
      stop('Backend warnings require review; no context SUCCESS')
    }
    path <- file.path(fitout,'all_results.tsv')
    if (!file.exists(path)) stop('Backend results missing')
    native <- read.delim(path,check.names=FALSE,stringsAsFactors=FALSE)
    required <- c('feature','metadata','coef','stderr','pval','qval')
    if (!all(required %in% names(native))) stop('Unexpected backend result schema')
    native <- native[native$metadata=='group',,drop=FALSE]
    if (anyDuplicated(native$feature) || !setequal(native$feature,ids[variable]))
      stop('Missing/duplicate/unexpected group results; review backend filtering/fits')
    if (any(!is.finite(native$coef)) || any(!is.finite(native$stderr)) ||
        any(native$stderr<=0) || any(!is.finite(native$pval)) ||
        any(native$pval<0 | native$pval>1)) stop('Invalid backend group inference')
    # MaAsLin2 may log convergence without issuing an R warning: inspect models.
    models_path <- file.path(fitout,'fits','models.rds')
    if (!file.exists(models_path)) stop('Saved models required for diagnostics')
    models <- readRDS(models_path)
    if (!is.list(models) || !length(models)) stop('Unexpected saved-model layout')
    for (model in models) {
      if (paired) {
        if (!inherits(model,'merMod')) stop('Expected paired mixed model')
        if (lme4::isSingular(model) || length(model@optinfo$conv$lme4$messages) ||
            any(model@optinfo$conv$opt!=0)) stop('Singular/nonconverged paired fit requires review')
      } else if (!inherits(model,'lm')) stop('Expected linear model')
    }
    i <- match(native$feature,ids)
    results$beta[i] <- native$coef
    results$stderr[i] <- native$stderr
    results$raw_p[i] <- native$pval
    results$native_q[i] <- native$qval
    results$estimable[i] <- TRUE
    results$status[i] <- 'ESTIMABLE'
  }
  bh <- full_family_bh(family,results$feature,results$raw_p,results$estimable)
  results$p_for_BH <- bh$p_for_BH
  results$wrapper_q <- bh$wrapper_q
  results$positive_discovery <- results$estimable & !is.na(results$beta) &
    results$beta>0 & results$wrapper_q<=.05
  write.table(results,file.path(output,'context_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  writeLines(c('status\tPASS_CONTEXT','validation\tDEVELOPMENT_NOT_DEFINITIVE'),file.path(output,'SUCCESS'))
  results
}
