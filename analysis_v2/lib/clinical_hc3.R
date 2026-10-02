# OLS coefficients unchanged; HC3 sandwich group standard error, t(n-p) reference.
# See https://sandwich.r-forge.r-project.org/reference/vcovHC.html
clinical_hc3 <- function(x,y) {
 if(!is.matrix(x)||!is.matrix(y)||nrow(x)!=nrow(y)||any(!is.finite(x))||any(!is.finite(y))||
    qr(x)$rank!=ncol(x)||nrow(x)<=ncol(x)||!'group' %in% colnames(x))stop('Invalid clinical design')
 j<-match('group',colnames(x));bread<-solve(crossprod(x));weights<-x%*%bread
 leverage<-rowSums(weights*x)
 if(any(1-leverage<1e-10))stop('Extreme leverage: HC3 undefined/unstable')
 fit<-lm.fit(x,y);beta<-matrix(fit$coefficients,nrow=ncol(x))[j,]
 r<-matrix(fit$residuals,nrow=nrow(x));df<-nrow(x)-ncol(x)
 ordinary_se<-sqrt(colSums(r^2)/df*bread[j,j])
 robust_se<-sqrt(colSums((r*(weights[,j]/(1-leverage)))^2))
 tolerance<-100*.Machine$double.eps*pmax(1,apply(abs(y),2,max))
 variable<-apply(y,2,function(v)diff(range(v)))>tolerance
 ok<-variable & apply(abs(r),2,max)>tolerance & is.finite(robust_se) & robust_se>0
 p<-classic<-rep(NA_real_,ncol(y))
 p[ok]<-2*pt(-abs(beta[ok]/robust_se[ok]),df)
 classic[ok]<-2*pt(-abs(beta[ok]/ordinary_se[ok]),df)
 data.frame(feature=colnames(y),beta=ifelse(ok,beta,NA_real_),estimable=ok,
   ordinary_se=ifelse(ok,ordinary_se,NA_real_),hc3_se=ifelse(ok,robust_se,NA_real_),
   ordinary_p=classic,hc3_p=p,ordinary_bh_q=p.adjust(ifelse(ok,classic,1),'BH'),
   hc3_bh_q=p.adjust(ifelse(ok,p,1),'BH'),residual_df=df,
   reference='HC3_WITH_RESIDUAL_DF_T_APPROXIMATION',stringsAsFactors=FALSE)
}
