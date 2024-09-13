library(EGG)
for(iters in 266:500){
    b=MASS::mvrnorm(m,rep(0,p),Sbb)*sqrt(n)/sqrt(m)
    u=MASS::mvrnorm(m,rep(0,p),Siguu)
    #if(UHP==1){
    #  b[1:round(m*UHPratio),1:5]=b[1:round(m*UHPratio),1:5]+UHPshift/sqrt(m)*sqrt(n)
    #}
    X=b+u
    remove(b,u)
    w=MASS::mvrnorm(10*m,rep(0,p),Siguu)
    Suu=cov(w)
    remove(w)
    fit1=EGG::entropy.mcp.spearman.sampling(BETA=X,Rnoise=Suu,lamvec=lamvec1*0.5,rho=0.05,subtime=300)
    ThetaList[iters,,]=fit1$Theta
    ThetaSEList[iters,,]=fit1$ThetaSE
    ThetaKList[iters,,]=fit1$K
 print(iters)
}


