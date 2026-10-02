args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: task repo')
task<-args[1];repo<-args[2]
source(file.path(repo,'analysis_v2/lib/clinical_hc3.R'))
source(file.path(repo,'analysis_v2/lib/clinical_wild_bootstrap.R'))
save<-function(x,name)write.table(x,file.path(task,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
s<-read.delim(file.path(task,'settings.tsv'),stringsAsFactors=FALSE)
a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE))
storage.mode(a)<-'double'
meta<-read.delim(file.path(task,'metadata.tsv'),row.names=1,stringsAsFactors=FALSE)
if(anyDuplicated(rownames(meta))||!setequal(rownames(meta),rownames(a))||any(!is.finite(a))||any(a<0))stop('Invalid input')
meta<-meta[rownames(a),,drop=FALSE]
if(anyDuplicated(meta$biological_sample_id))stop('Independent people required')
meta$sex<-factor(meta$sex,levels=c('Female','Male'));meta$age<-as.numeric(meta$age)
if(anyNA(meta[,c('group','age','sex')])||!setequal(meta$group,c(0,1)))stop('Invalid clinical metadata')
x<-model.matrix(~group+age+sex,meta)
# Lexical selection is fixed independently of outcomes; subset q-values are not full-family q-values.
features<-head(sort(colnames(a)),s$family_n)
if(length(features)!=s$family_n)stop('Insufficient engineering features')
save(data.frame(feature=features,selection='FIRST_20_LEXICAL_NOT_OUTCOME_SELECTED'),'engineering_family.tsv')
if(s$scenario=='actual'){
 y<-log2(1+a[,features,drop=FALSE]/1e-8)
 fit<-clinical_wild_bootstrap(x,y,s$draws,20261003L)
 fit$interpretation<-'ENGINEERING_SUBSET_Q_NOT_FULL_FAMILY_Q'
 save(fit,'comparison.tsv')
}else{
 draws<-lapply(seq_len(s$simulations),function(k){
  i<-s$chunk*s$simulations+k;set.seed(20261003L+i)
  n<-nrow(x);m<-length(features)
  e<-switch(s$scenario,gaussian=matrix(rnorm(n*m),n),skewed=matrix(rexp(n*m)-1,n),heteroskedastic=matrix(rnorm(n*m),n)*ifelse(meta$group==1,3,1),stop('Unknown scenario'))
  y<-e+.3*as.numeric(scale(meta$age))+.2*as.numeric(meta$sex=='Male')
  dimnames(y)<-list(rownames(x),features)
  f<-clinical_wild_bootstrap(x,y,s$draws,300000L+i)
  data.frame(simulation=i,method=c('ordinary','HC3','wild'),discoveries=c(sum(f$ordinary_bh_q<=.05),sum(f$hc3_bh_q<=.05),sum(f$wild_bh_q<=.05)),pointwise_rejections=c(sum(f$ordinary_p<=.05),sum(f$hc3_p<=.05),sum(f$wild_p<=.05)),family_n=m)
 });save(do.call(rbind,draws),'draws.tsv')
}
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
