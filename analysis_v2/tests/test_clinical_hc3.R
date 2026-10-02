source('analysis_v2/lib/clinical_hc3.R')
set.seed(19);meta<-data.frame(group=rep(0:1,each=20),age=rnorm(40),sex=factor(rep(c('Female','Male'),20)))
x<-model.matrix(~group+age+sex,meta);y<-cbind(a=rnorm(40)*ifelse(meta$group==1,3,1),b=rnorm(40),constant=1)
f<-clinical_hc3(x,y)
for(j in 1:2){
 fit<-lm(y[,j]~group+age+sex,data=meta);s<-summary(fit)$coefficients
 stopifnot(isTRUE(all.equal(f$beta[j],unname(s['group','Estimate']),tolerance=1e-10)))
 bread<-solve(crossprod(x));e<-residuals(fit);h<-hatvalues(fit)
 v<-bread%*%crossprod(x,x*(e/(1-h))^2)%*%bread
 stopifnot(isTRUE(all.equal(f$hc3_se[j],sqrt(v['group','group']),tolerance=1e-10)))
 if(requireNamespace('sandwich',quietly=TRUE))stopifnot(isTRUE(all.equal(f$hc3_se[j],sqrt(sandwich::vcovHC(fit,type='HC3')['group','group']),tolerance=1e-10)))
}
stopifnot(!f$estimable[3],is.na(f$hc3_p[3]),f$hc3_bh_q[3]==1)
if(!requireNamespace('sandwich',quietly=TRUE))cat('[NOTE] Package comparison unavailable; independent matrix formula passed\n')
cat('[PASS] HC3 matrix formula, OLS coefficient equivalence, degeneracy\n')
