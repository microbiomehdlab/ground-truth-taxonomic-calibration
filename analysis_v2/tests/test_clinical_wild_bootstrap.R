source('analysis_v2/lib/clinical_hc3.R')
source('analysis_v2/lib/clinical_wild_bootstrap.R')
set.seed(42)
n<-40L
x<-cbind(intercept=1,group=rep(0:1,each=n/2),age=rnorm(n),sex=rep(0:1,n/2))
rownames(x)<-paste0('person_',seq_len(n))
y<-cbind(null=rnorm(n),signal=8*x[,'group']+rnorm(n),constant=1)
rownames(y)<-rownames(x)
rng<-.Random.seed
a<-clinical_wild_bootstrap(x,y,199L,17L)
stopifnot(identical(rng,.Random.seed),a$wild_p[2]<=.01,
          is.na(a$wild_p[3]),a$wild_bh_q[3]==1,
          identical(a,clinical_wild_bootstrap(x,y,199L,17L)))
shard<-clinical_wild_bootstrap(x,y[,2,drop=FALSE],199L,17L)
stopifnot(a$wild_p[2]==shard$wild_p, a$beta[2]==shard$beta)
shifted<-y+3*x[,'age']-2*x[,'sex']
b<-clinical_wild_bootstrap(x,shifted,199L,17L)
stopifnot(all.equal(a$wild_p[1:2],b$wild_p[1:2]),
          isTRUE(all.equal(a$beta[1:2],b$beta[1:2])))
stopifnot(inherits(try(clinical_wild_bootstrap(x,y[n:1,,drop=FALSE],19L),silent=TRUE),'try-error'))
# Independent scalar lm/vcov sandwich calculation, explicit weights and restricted fit.
weights<-matrix(sample(c(-1,1),n*39,replace=TRUE),n)
one<-y[,1,drop=FALSE]
tested<-clinical_wild_bootstrap(x,one,39L,multipliers=weights)
x0<-x[,-2,drop=FALSE]; restricted<-lm(one[,1]~x0-1)
h0<-hatvalues(restricted); e0<-residuals(restricted)/sqrt(1-h0)
statistic<-function(z){
 f<-lm(z~x-1);inv<-solve(crossprod(x));h<-hatvalues(f)
 covariance<-inv%*%crossprod(x,x*(residuals(f)/(1-h))^2)%*%inv
 abs(coef(f)[2]/sqrt(covariance[2,2]))
}
threshold<-statistic(one[,1])
oracle<-(1+sum(vapply(seq_len(39),function(k)
 statistic(fitted(restricted)+e0*weights[,k])>=threshold,logical(1))))/40
stopifnot(tested$wild_p==oracle)
cat('[PASS] Determinism, sample identities, nuisance shift, signal, constants, feature sharding and RNG preservation\n')
