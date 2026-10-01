# Pure input/transformation/multiplicity contracts; no package or model run.
prepare_maaslin_input <- function(abundance, metadata, family, paired=FALSE,
                                  pseudocount=1e-8) {
  if (!is.matrix(abundance) || !is.numeric(abundance)) stop('Numeric matrix required')
  if (is.null(rownames(abundance)) || is.null(colnames(abundance)) ||
      anyDuplicated(rownames(abundance)) || anyDuplicated(colnames(abundance)))
    stop('Unique sample and feature identifiers required')
  if (!nrow(abundance) || !ncol(abundance) ||
      any(!is.finite(abundance)) || any(abundance < 0 | abundance > 1.00001))
    stop('Finite native abundance fractions required')
  if (!is.data.frame(metadata) || is.null(rownames(metadata)) ||
      anyDuplicated(rownames(metadata)) ||
      !setequal(rownames(abundance), rownames(metadata)))
    stop('Exact feature/metadata sample identity required; no silent intersection')
  if (!length(family) || anyNA(family) || any(!nzchar(family)) || anyDuplicated(family) ||
      !all(family %in% colnames(abundance))) stop('Frozen family coverage required')
  if (length(pseudocount)!=1 || !is.finite(pseudocount) || pseudocount<=0)
    stop('Positive finite pseudocount required')
  metadata <- metadata[rownames(abundance),,drop=FALSE]
  if (paired) {
    if (!all(c('biological_sample_id','spike_state') %in% names(metadata)))
      stop('Pairing identity and spike_state required')
    if (anyNA(metadata$biological_sample_id) || any(!nzchar(as.character(metadata$biological_sample_id))) ||
        anyNA(metadata$spike_state) || !all(metadata$spike_state %in% c('original','spiked')))
      stop('Invalid pairing metadata')
    pairs <- table(metadata$biological_sample_id,
                   factor(metadata$spike_state,levels=c('original','spiked')))
    if (any(pairs!=1)) stop('Exactly one original and one spiked observation per person required')
  }
  native <- abundance[,family,drop=FALSE]
  # Algebraically log2(a+p)-log2(p), with zeros exactly zero.
  transformed <- log2(1 + native/pseudocount)
  if (any(!is.finite(transformed))) stop('Transformation overflow')
  list(data=as.data.frame(transformed,check.names=FALSE), metadata=metadata,
       native=native, family=family, pseudocount=pseudocount,
       normalization='NONE', transform='NONE', standardize=FALSE)
}

full_family_bh <- function(family, feature, pvalue, estimable) {
  if (!length(family) || anyNA(family) || anyDuplicated(family) ||
      any(!nzchar(family))) stop('Invalid family')
  if (length(feature)!=length(pvalue) || length(feature)!=length(estimable) ||
      anyNA(feature) || anyDuplicated(feature) || !setequal(family,feature))
    stop('Results must account for every frozen family member exactly once')
  if (!is.logical(estimable) || anyNA(estimable)) stop('Explicit estimability required')
  if (any(!is.finite(pvalue[estimable])) || any(pvalue[estimable]<0 | pvalue[estimable]>1))
    stop('Invalid estimable p value')
  if (any(!is.na(pvalue[!estimable]))) stop('Non-estimable raw p must be NA')
  order <- match(family,feature)
  p <- pvalue[order]
  ok <- estimable[order]
  bookkeeping <- ifelse(ok,p,1)
  data.frame(feature=family,raw_p=p,estimable=ok,p_for_BH=bookkeeping,
             wrapper_q=p.adjust(bookkeeping,method='BH'),stringsAsFactors=FALSE)
}
