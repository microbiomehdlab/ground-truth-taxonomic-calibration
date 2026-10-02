#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop('Usage: task-directory repository')
task <- args[1]; repo <- args[2]
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/paired_difference_context.R'))
context <- read.delim(file.path(task,'context.tsv'),check.names=FALSE)
if (nrow(context)!=1 || context$analysis!='DA2' || context$paired!=1)
  stop('Paired-difference inference is restricted to paired DA2')
data <- as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE))
storage.mode(data) <- 'double'
meta <- read.delim(file.path(task,'metadata.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
answer <- paired_difference_results(data,meta,colnames(data))
out <- file.path(task,'fit')
if (file.exists(out)) stop('Fresh output required')
dir.create(out)
write.table(answer$results,file.path(out,'context_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(data.frame(biological_sample_id=rownames(answer$differences),answer$differences,check.names=FALSE),
  file.path(out,'paired_differences.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c('method\tpaired_difference_t','pseudocount\t1e-8','alternative\ttwo.sided',
  'multiplicity\tfrozen_full_family_BH','definitive\t0'),file.path(out,'method.tsv'))
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
writeLines('status\tPASS_ENGINEERING_CONTEXT',file.path(out,'SUCCESS'))
