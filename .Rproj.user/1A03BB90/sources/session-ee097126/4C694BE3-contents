source("basicfunction.R")
source("mainfunction.R")
for(iters in iterstar:Miter){
results <- foreach(i=c(1:(numCores)),.errorhandling="pass",.combine=c) %dopar% {
library(clime)
b=MASS::mvrnorm(m,rep(0,p),Sbb)*sqrt(n)/sqrt(m)
u=MASS::mvrnorm(m,rep(0,p),Siguu)
if(UHP==1){
b[1:round(m*UHPratio),1:5]=b[1:round(m*UHPratio),1:5]+UHPshift/sqrt(m)*sqrt(n)
}
X=b+u
remove(b,u)
w=MASS::mvrnorm(10*m,rep(0,p),Siguu)
Suu=cov(w)
remove(w)
fit11=entropy.mcp(X,Suu,lamvec=lamvec1,rho=0.05,method="spearman",IC="BIC",thres=0.75)
fit1=entropy.mcp(X,Suu,lamvec=lamvec1,rho=0.05,method="pearson",IC="BIC",thres=0.75,Theta.ini=fit11$Theta)
fit22=graph.lasso(X,lamvec=lamvec1,method="spearman",IC="BIC",thres=0.75)
fit2=graph.lasso(X,lamvec=lamvec1,method="pearson",IC="BIC",thres=0.75)
fit33=dtrace.lasso(X,Suu*0,lamvec=lamvec1*1.5,rho=0.05,method="spearman",IC="BIC",thres=0.75,Theta.ini=fit11$Theta)
fit3=dtrace.lasso(X,Suu*0,lamvec=lamvec1*1.5,rho=0.05,method="pearson",IC="BIC",thres=0.75,Theta.ini=fit11$Theta)
fit4=clime(X,nlambda=50,perturb=F,lambda.max=0.5,lambda.min=0.01)
fit4=climeselect(X,fit4)
fit44=clime(normtransform(X),nlambda=50,perturb=F,lambda.max=0.5,lambda.min=0.01)
fit44=climeselect(normtransform(X),fit44)

error1=c(entropyloss(Rbb,fit1$Theta),quadraticloss(Rbb,fit1$Theta))
error2=c(entropyloss(Rbb,fit2$Theta),quadraticloss(Rbb,fit2$Theta))
error3=c(entropyloss(Rbb,fit3$Theta),quadraticloss(Rbb,fit3$Theta))
error4=c(entropyloss(Rbb,fit4$Theta),quadraticloss(Rbb,fit4$Theta))
error11=c(entropyloss(Rbb,fit11$Theta),quadraticloss(Rbb,fit11$Theta))
error22=c(entropyloss(Rbb,fit22$Theta),quadraticloss(Rbb,fit22$Theta))
error33=c(entropyloss(Rbb,fit33$Theta),quadraticloss(Rbb,fit33$Theta))
error44=c(entropyloss(Rbb,fit44$Theta),quadraticloss(Rbb,fit44$Theta))
tptn1=TPTN(Tbb,fit1$Theta)
tptn2=TPTN(Tbb,fit2$Theta)
tptn3=TPTN(Tbb,fit3$Theta)
tptn4=TPTN(Tbb,fit4$Theta)
tptn11=TPTN(Tbb,fit11$Theta)
tptn22=TPTN(Tbb,fit22$Theta)
tptn33=TPTN(Tbb,fit33$Theta)
tptn44=TPTN(Tbb,fit44$Theta)
return(c(error1,tptn1,error2,tptn2,error3,tptn3,error4,tptn4,error11,tptn11,error22,tptn22,error33,tptn33,error44,tptn44))
}
B=array(results,c(4,8,numCores))
ind=c(((iters-1)*numCores+1):(iters*numCores))
RES1[ind,]=t(B[1,,])
RES2[ind,]=t(B[2,,])
RES3[ind,]=t(B[3,,])
RES4[ind,]=t(B[4,,])

print(iters)
print(colMeans(RES1[1:max(ind),]))
print(colMeans(RES2[1:max(ind),]))
print(colMeans(RES3[1:max(ind),]))
print(colMeans(RES4[1:max(ind),]))
save.image(filenames)
}



