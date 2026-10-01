#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop('Usage: run_da_pilot_context.R task-directory repository')
task <- args[1]; repo <- args[2]
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/maaslin_context.R'))
context <- read.delim(file.path(task,'context.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
if (nrow(context)!=1) stop('One context required')
data <- as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE))
storage.mode(data) <- 'double'
meta <- read.delim(file.path(task,'metadata.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
meta$group <- as.integer(meta$group)
covariates <- if (context$analysis=='DA1') c('age','sex') else character()
if (length(covariates)) {
    meta$age <- as.numeric(meta$age)
    meta$sex <- factor(meta$sex,levels=c('Female','Male'))
}
if (!requireNamespace('Maaslin2',quietly=TRUE) || as.character(packageVersion('Maaslin2'))!='1.18.0')
    stop('Require pinned MaAsLin2 1.18.0')
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
fit_maaslin_context(data,meta,colnames(data),file.path(task,'fit'),
                   paired=context$paired==1,covariates=covariates)
