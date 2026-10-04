args <- commandArgs(TRUE); stopifnot(length(args)==4)
task <- args[1]; scratch <- args[2]; out <- args[3]; repo <- args[4]
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/paired_difference_context.R'))
contexts <- jsonlite::fromJSON(file.path(task,'contexts.json'),simplifyVector=FALSE)
targets <- jsonlite::fromJSON(file.path(task,'targets.json'))
full <- gzfile(file.path(out,'full_family_results.tsv.gz'),'wt')
first <- TRUE; selected <- list()
for (context in contexts) {
  id <- context$context_id; folder <- file.path(scratch,id)
  x <- as.matrix(read.delim(file.path(folder,'abundance.tsv'),row.names=1,check.names=FALSE)); storage.mode(x)<-'double'
  m <- read.delim(file.path(folder,'metadata.tsv'),row.names=1,check.names=FALSE)
  family <- unlist(context$family)
  result <- paired_difference_results(x,m,family)$results
  result$significant_two_sided <- result$estimable & result$wrapper_q<=.05
  result <- cbind(data.frame(context_id=id,cohort=context$cohort,profiler=context$profiler,
      condition=context$condition,population=context$population,spike_label=context$spike_label,
      nominal_total_dose=context$nominal_total_dose,family_n=length(family)),result)
  write.table(result,full,sep='\t',row.names=FALSE,col.names=first,quote=FALSE,na='NA'); first<-FALSE
  mapping <- unlist(targets[[context$profiler]])
  keep <- result$feature %in% names(mapping)
  if (context$population=='independent') keep <- keep & result$feature %in% names(mapping)[mapping==context$spike_label]
  rr <- result[keep,,drop=FALSE]; rr$target_label<-unname(unlist(mapping[rr$feature]))
  selected[[length(selected)+1]] <- rr
  cat('[PASS] Paired supplementary context:',id,'\n')
}
close(full)
write.table(do.call(rbind,selected),file.path(out,'target_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
