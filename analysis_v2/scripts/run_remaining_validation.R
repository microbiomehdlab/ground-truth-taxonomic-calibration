args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3)stop('Usage: task mode repo')
task<-args[1];mode<-args[2];repo<-args[3]
source(file.path(repo,'analysis_v2/lib/maaslin_contract.R'))
source(file.path(repo,'analysis_v2/lib/maaslin_context.R'))
source(file.path(repo,'analysis_v2/lib/exact_da3_candidate.R'))
source(file.path(repo,'analysis_v2/lib/monte_carlo_da3_candidate.R'))
save<-function(x,name)write.table(x,file.path(task,name),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
settings<-read.delim(file.path(task,'settings.tsv'),stringsAsFactors=FALSE)
if(mode %in% c('clinical_backend','clinical_null')){
 a<-as.matrix(read.delim(file.path(task,'abundance.tsv'),row.names=1,check.names=FALSE));storage.mode(a)<-'double'
 meta<-read.delim(file.path(task,'metadata.tsv'),row.names=1,stringsAsFactors=FALSE)
 if(!setequal(rownames(a),rownames(meta)))stop('Sample mismatch')
 meta<-meta[rownames(a),,drop=FALSE];meta$age<-as.numeric(meta$age);meta$sex<-factor(meta$sex,levels=c('Female','Male'))
 if(anyNA(meta[,c('group','age','sex')]))stop('Invalid clinical covariates')
 x<-model.matrix(~group+age+sex,meta)
 if(qr(x)$rank!=ncol(x)||nrow(x)<=ncol(x))stop('Unidentifiable clinical design')
 calculate<-function(y){
   f<-lm.fit(x,y);df<-nrow(x)-ncol(x);v<-colSums(f$residuals^2)/df
   se<-sqrt(v*solve(crossprod(x))[2,2]);beta<-f$coefficients[2,]
   p<-2*pt(-abs(beta/se),df);list(beta=beta,p=p)
 }
 if(mode=='clinical_backend'){
   fit<-fit_maaslin_context(a,meta,colnames(a),file.path(task,'backend'),covariates=c('age','sex'))
   independent<-calculate(log2(1+a/1e-8));ok<-fit$estimable
   if(!isTRUE(all.equal(fit$beta[ok],unname(independent$beta[ok]),tolerance=1e-7))||
      !isTRUE(all.equal(fit$raw_p[ok],unname(independent$p[ok]),tolerance=1e-7)))stop('Clinical backend disagreement')
   save(data.frame(checked_features=sum(ok),design_rank=ncol(x),residual_df=nrow(x)-ncol(x),status='PASS_PINNED_CLINICAL_BACKEND'),'summary.tsv')
 }else{
   set.seed(20261003);m<-ncol(a);n<-nrow(a);scenario<-settings$scenario
   draws<-lapply(seq_len(200),function(i){
    e<-switch(scenario,gaussian=matrix(rnorm(n*m),n),skewed=matrix(rexp(n*m)-1,n),heteroskedastic=matrix(rnorm(n*m),n)*ifelse(meta$group==1,3,1),stop('Unknown scenario'))
    y<-e+as.numeric(scale(meta$age))*.3+as.numeric(meta$sex=='Male')*.2
    fit<-calculate(y);data.frame(repetition=i,discoveries=sum(p.adjust(fit$p,'BH')<=.05))
   });draws<-do.call(rbind,draws);save(draws,'draws.tsv')
   ci<-binom.test(sum(draws$discoveries>0),nrow(draws))$conf.int
   save(data.frame(scenario=scenario,repetitions=200,any_discovery_rate=mean(draws$discoveries>0),ci_low=ci[1],ci_high=ci[2],interpretation='SYNTHETIC_CONDITIONAL_ON_ACTUAL_DESIGN_NOT_CLINICAL_ERROR_PROOF'),'summary.tsv')
 }
}else if(mode=='partial_null'){
 n<-settings$n;m<-settings$family_n;scenario<-settings$scenario;g<-rep(0:1,each=n);set.seed(20261003+n+m)
 draws<-lapply(seq_len(200),function(i){
  e<-switch(scenario,gaussian=matrix(rnorm(2*n*m),2*n),correlated=matrix(rnorm(2*n*m),2*n)*sqrt(.7)+rnorm(2*n)*sqrt(.3),skewed=matrix(rexp(2*n*m)-1,2*n),stop('Unknown scenario'))
  e[g==1,1:10]<-e[g==1,1:10]+1
  colnames(e)<-sprintf('f%04d',seq_len(m))
  fit<-if(n==5)exact_da3_candidate(e,g)else monte_carlo_da3_candidate(e,g,repetitions=9999,seed=20261003+i)
  # Only features 11..m have a known unchanged distribution in this synthetic design.
  data.frame(repetition=i,false_discoveries=sum(fit$max_statistic_fwer_p[-(1:10)]<=.05),true_positive_discoveries=sum(fit$beta[1:10]>0 & fit$max_statistic_fwer_p[1:10]<=.05))
 });draws<-do.call(rbind,draws);save(draws,'draws.tsv');ci<-binom.test(sum(draws$false_discoveries>0),200)$conf.int
 save(data.frame(n=n,family_n=m,scenario=scenario,repetitions=200,false_any_discovery_rate=mean(draws$false_discoveries>0),ci_low=ci[1],ci_high=ci[2],mean_true_positives=mean(draws$true_positive_discoveries),interpretation='PARTIAL_NULL_SYNTHETIC_NOT_ARBITRARY_STRONG_FWER_PROOF'),'summary.tsv')
}else stop('Unknown mode')
writeLines(capture.output(sessionInfo()),file.path(task,'sessionInfo.txt'))
