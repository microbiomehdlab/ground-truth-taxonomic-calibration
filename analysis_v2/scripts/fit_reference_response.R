# Direct pinned MaAsLin2 reference fit; observed primary results are reused.
fit_reference_response <- function(input, output, repo, backend=NULL) {
  for (name in c('maaslin_contract.R','maaslin_context.R','clinical_hc3.R'))
    source(file.path(repo,'analysis_v2/lib',name))
  if (is.null(backend) && (!requireNamespace('Maaslin2',quietly=TRUE) ||
      as.character(packageVersion('Maaslin2'))!='1.18.0')) stop('Require MaAsLin2 1.18.0')
  if (dir.exists(output)) stop('Fresh output required')
  contexts <- read.delim(file.path(input,'contexts.tsv'),stringsAsFactors=FALSE)
  if (!nrow(contexts) || anyDuplicated(contexts$context_id)) stop('Invalid contexts')
  dir.create(output,recursive=TRUE)
  full <- gzfile(file.path(output,'features.tsv.gz'),'wt'); on.exit(close(full),add=TRUE)
  target_rows <- list(); first <- TRUE
  for (id in contexts$context_id) {
    folder <- file.path(input,id)
    metadata <- read.delim(file.path(folder,'metadata.tsv'),row.names=1,check.names=FALSE)
    targets <- read.delim(file.path(folder,'targets.tsv'),stringsAsFactors=FALSE)
    for (arm in c('reference','observed')) {
      a <- as.matrix(read.delim(file.path(folder,paste0(arm,'.tsv')),row.names=1,check.names=FALSE))
      storage.mode(a) <- 'double'
      if (!identical(rownames(a),rownames(metadata)) || anyDuplicated(rownames(a)) ||
          any(!is.finite(a)) || any(a<0) || any(a>1.00001) ||
          anyDuplicated(targets$feature) || !all(targets$feature %in% colnames(a))) stop('Input identity/value mismatch')
      n <- contexts$n[match(id,contexts$context_id)]
      if(nrow(a)!=2*n || any(table(factor(metadata$group,levels=0:1))!=n)) stop('Group size differs')
      y <- log2(1+a/1e-8)
      hc <- clinical_hc3(model.matrix(~group,metadata),y)
      if (arm=='reference') {
        fit <- fit_maaslin_context(a,metadata,colnames(a),file.path(folder,'reference_fit'),backend=backend)
      } else {
        # No replacement observed primary fit: this arm only supplies sensitivity.
        fit <- data.frame(feature=colnames(a),beta=hc$beta,stderr=hc$ordinary_se,
          raw_p=hc$ordinary_p,wrapper_q=hc$ordinary_bh_q,estimable=hc$estimable,
          status=ifelse(hc$estimable,'ESTIMABLE','NON_ESTIMABLE_HC3'),
          positive_discovery=hc$estimable & !is.na(hc$beta) & hc$beta>0 & hc$ordinary_bh_q<=.05)
      }
      if(!identical(fit$feature,hc$feature)) stop('Full family order differs')
      fit$hc3_p <- hc$hc3_p; fit$hc3_q <- hc$hc3_bh_q
      fit$hc3_stderr <- hc$hc3_se; fit$hc3_estimable <- hc$estimable
      fit$context_id <- id; fit$arm <- arm
      keep <- c('context_id','arm','feature','beta','stderr','raw_p','wrapper_q','status',
                'estimable','positive_discovery','hc3_p','hc3_q','hc3_stderr','hc3_estimable')
      fit <- fit[,keep]
      write.table(fit,full,sep='\t',quote=FALSE,row.names=FALSE,col.names=first,na='NA'); first <- FALSE
      selected <- fit[match(targets$feature,fit$feature),,drop=FALSE]
      selected$target_label <- targets$target_label
      target_rows[[length(target_rows)+1L]] <- selected
    }
    cat('[PASS] Matched reference/sensitivity:',id,'\n')
  }
  write.table(do.call(rbind,target_rows),file.path(output,'targets.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))
}
if(sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if(length(args)!=3) stop('Usage: INPUT OUTPUT REPO')
  fit_reference_response(args[1],args[2],args[3])
}
