# Same direct backend as the timing trial; compact retention, no new model.
fit_maaslin_batch <- function(input, out, repo, backend=NULL) {
  for (f in c('maaslin_contract.R','maaslin_context.R','unpaired_null.R'))
    source(file.path(repo,'analysis_v2/lib',f))
  source(file.path(repo,'analysis_v2/scripts/time_maaslin_context.R'))
  if (is.null(backend)) {
    if (!requireNamespace('Maaslin2',quietly=TRUE) ||
        as.character(packageVersion('Maaslin2'))!='1.18.0') stop('Require MaAsLin2 1.18.0')
  }
  tasks <- read.delim(file.path(input,'contexts.tsv'),check.names=FALSE,stringsAsFactors=FALSE,
                     colClasses='character',na.strings=character())
  if (!nrow(tasks) || anyDuplicated(tasks$context_id)) stop('Invalid contexts')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  connections <- lapply(c('features.tsv.gz','targets.tsv.gz','differences.tsv.gz'),
                        function(n) gzfile(file.path(out,n),'wt'))
  on.exit(lapply(connections,close),add=TRUE)
  summaries <- list(); first_difference <- TRUE
  for (i in seq_len(nrow(tasks))) {
    context <- tasks[i,,drop=FALSE]; context$n <- as.integer(context$n); id <- context$context_id
    folder <- file.path(input,id)
    writeLines(id,file.path(out,'current_context.txt'))
    a <- as.matrix(read.delim(file.path(folder,'abundance.tsv'),row.names=1,check.names=FALSE))
    storage.mode(a) <- 'double'
    m <- read.delim(file.path(folder,'metadata.tsv'),row.names=1,check.names=FALSE)
    targets <- read.delim(file.path(folder,'targets.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
    if (!identical(rownames(a),rownames(m)) || nrow(a)!=2*context$n ||
        any(table(factor(m$group,levels=0:1))!=context$n) || anyDuplicated(rownames(a)))
      stop('Independent-group identity/count mismatch')
    if (anyDuplicated(targets$feature) || anyDuplicated(targets$target_label) ||
        !all(targets$feature %in% colnames(a))) stop('Invalid targets')
    start <- proc.time()[['elapsed']]
    fit <- fit_maaslin_context(a,m,colnames(a),file.path(folder,'native_fit'),backend=backend)
    fast <- unpaired_null_results(log2(1+a/1e-8),m$group)
    comparison <- compare_maaslin_timing(fit,fast)
    if (!identical(fit$feature,colnames(a))) stop('Full-family output order mismatch')
    fit <- data.frame(context_id=id,fit,stringsAsFactors=FALSE)
    write.table(fit,connections[[1]],sep='\t',quote=FALSE,row.names=FALSE,col.names=i==1,na='NA')
    selected <- fit[match(targets$feature,fit$feature),,drop=FALSE]
    selected$target_label <- targets$target_label
    write.table(selected,connections[[2]],sep='\t',quote=FALSE,row.names=FALSE,col.names=i==1,na='NA')
    differences <- comparison[!comparison$agreement,,drop=FALSE]
    if (nrow(differences) || first_difference) {
      write.table(data.frame(context_id=rep(id,nrow(differences)),differences),connections[[3]],sep='\t',quote=FALSE,
                  row.names=FALSE,col.names=first_difference,na='NA')
      first_difference <- FALSE
    }
    summaries[[i]] <- data.frame(context, family_n=nrow(fit), estimable=sum(fit$estimable),
        positive_discoveries=sum(fit$positive_discovery), target_rows=nrow(selected),
        mismatched_features=nrow(differences), seconds=proc.time()[['elapsed']]-start,
        backend_version=if(is.null(backend)) as.character(packageVersion('Maaslin2')) else 'TEST_DOUBLE')
    # Only this newly created scratch context is removed, after diagnostic checks and writes.
    unlink(file.path(folder,'native_fit'),recursive=TRUE)
    cat('[PASS] Direct MaAsLin2:',id,'features:',nrow(fit),'\n')
  }
  write.table(do.call(rbind,summaries),file.path(out,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (length(args)!=3) stop('Usage: fit_maaslin_batch.R INPUT OUTPUT REPO')
  fit_maaslin_batch(args[1],args[2],args[3])
}
