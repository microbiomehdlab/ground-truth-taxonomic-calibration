args <- commandArgs(TRUE)
repo <- if(length(args)) args[1] else getwd()
source(file.path(repo,'analysis_v2/scripts/fit_reference_response.R'))
temporary <- tempfile('reference-fixture-'); dir.create(temporary)
input <- file.path(temporary,'input');dir.create(input)
id <- 'da3_0000000'; folder <- file.path(input,id); dir.create(folder)
ids <- paste0('sample',1:20)
group <- rep(0:1,each=10)
a <- matrix(.001*(1:20),nrow=20,ncol=12)
colnames(a) <- c(paste0('taxon',1:10),'constant','perfect');rownames(a)<-ids
a[,'constant'] <- 0;a[,'perfect']<-group*.1
reference <- a;reference[,1:10] <- reference[,1:10]+group*.01
write.table(data.frame(observation_id=ids,group=group),file.path(folder,'metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(arm in c('reference','observed'))write.table(data.frame(observation_id=ids,if(arm=='reference')reference else a,check.names=FALSE),file.path(folder,paste0(arm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(feature=paste0('taxon',1:10),target_label=paste0('L',1:10)),file.path(folder,'targets.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(context_id=id,n=10),file.path(input,'contexts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Test double checks orchestration; --actual exercises the pinned real backend.
backend <- function(input_data,input_metadata,output,...) {
  dir.create(file.path(output,'fits'),recursive=TRUE)
  models <- lapply(input_data,function(y)lm(y~group,data=input_metadata))
  rows <- lapply(seq_along(models),function(i) {
    co <- summary(models[[i]])$coefficients['group',]
    data.frame(feature=names(input_data)[i],metadata='group',coef=co[1],stderr=co[2],pval=co[4],qval=1)
  })
  write.table(do.call(rbind,rows),file.path(output,'all_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  saveRDS(models,file.path(output,'fits/models.rds'))
}
fit_reference_response(input,file.path(temporary,'output'),repo,backend=if('--actual'%in%args)NULL else backend)
r <- read.delim(file.path(temporary,'output/targets.tsv'))
f <- read.delim(gzfile(file.path(temporary,'output/features.tsv.gz')))
stopifnot(nrow(r)==20,nrow(f)==24,all(table(r$arm)==10),all(f$wrapper_q>=0 & f$wrapper_q<=1),
          all(!f$estimable[f$feature=='constant']),all(!f$estimable[f$feature=='perfect']),
          all(f$wrapper_q[f$feature=='constant']==1))
cat('[PASS] Full family, direct reference backend, HC3 on both arms and non-estimable retention\n')
