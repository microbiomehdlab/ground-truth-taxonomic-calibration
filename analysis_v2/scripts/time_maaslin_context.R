# Direct pinned MaAsLin2 versus the existing DA3 group-only LM calculation.
compare_maaslin_timing <- function(actual, fast) {
  if (anyDuplicated(actual$feature) || anyDuplicated(fast$feature) ||
      !setequal(actual$feature, fast$feature)) stop('Feature identity mismatch')
  fast <- fast[match(actual$feature, fast$feature), , drop=FALSE]
  near <- function(a,b) {
    same_na <- is.na(a) & is.na(b)
    same_na | (!is.na(a) & !is.na(b) & is.finite(a) & is.finite(b) &
                 abs(a-b) <= 1e-8 + 1e-6*pmax(abs(a),abs(b)))
  }
  positive <- fast$estimable & !is.na(fast$beta) & fast$beta>0 & fast$wrapper_q<=.05
  out <- data.frame(feature=actual$feature, actual_status=actual$status, fast_status=fast$status,
    actual_beta=actual$beta, fast_beta=fast$beta, actual_p=actual$raw_p, fast_p=fast$raw_p,
    native_q=actual$native_q, actual_full_family_q=actual$wrapper_q, fast_full_family_q=fast$wrapper_q,
    status_equal=actual$status==fast$status & actual$estimable==fast$estimable,
    beta_equal=near(actual$beta,fast$beta), p_equal=near(actual$raw_p,fast$raw_p),
    q_equal=near(actual$wrapper_q,fast$wrapper_q),
    discovery_equal=actual$positive_discovery==positive, stringsAsFactors=FALSE)
  out$agreement <- with(out,status_equal & beta_equal & p_equal & q_equal & discovery_equal)
  if (anyNA(out$agreement)) stop('Undefined comparison')
  out
}

time_maaslin_context <- function(folder, repo, backend=NULL) {
  for (name in c('maaslin_contract.R','maaslin_context.R','unpaired_null.R'))
    source(file.path(repo,'analysis_v2','lib',name))
  a <- read.delim(file.path(folder,'abundance.tsv'), check.names=FALSE, row.names=1)
  m <- read.delim(file.path(folder,'metadata.tsv'), check.names=FALSE, row.names=1)
  a <- as.matrix(a); storage.mode(a) <- 'double'
  input <- prepare_maaslin_input(a,m,colnames(a),paired=FALSE)
  m <- input$metadata
  if (any(table(factor(m$group,levels=0:1))<2)) stop('Too few independent people')
  repetitions <- 50L
  start <- proc.time()
  for (i in seq_len(repetitions)) fast <- unpaired_null_results(log2(1+a/1e-8),m$group)
  fast_time <- (proc.time()-start)/repetitions
  start <- proc.time()
  actual <- fit_maaslin_context(a,m,colnames(a),file.path(folder,'direct_maaslin'),backend=backend)
  native_time <- proc.time()-start
  comparison <- compare_maaslin_timing(actual,fast)
  write.table(comparison,file.path(folder,'comparison.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  timing <- data.frame(maaslin_version=if(is.null(backend)) as.character(packageVersion('Maaslin2')) else 'TEST_DOUBLE',
    R_version=as.character(getRversion()), maaslin_seconds=unname(native_time['elapsed']),
    maaslin_cpu_seconds=unname(sum(native_time[c('user.self','sys.self')])),
    fast_seconds=unname(fast_time['elapsed']), fast_repetitions=repetitions,
    backend_features=sum(actual$estimable), mismatched_features=sum(!comparison$agreement),
    discovery_disagreements=sum(!comparison$discovery_equal))
  write.table(timing,file.path(folder,'timing.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(folder,'sessionInfo.txt'))
  cat('[PASS] Direct MaAsLin2 computation; mismatched features:',sum(!comparison$agreement),'\n')
  invisible(comparison)
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (length(args)!=2) stop('Usage: time_maaslin_context.R FOLDER REPO')
  time_maaslin_context(args[1],args[2])
}
