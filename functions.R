source("/basicfunction.R")

entropyloss=function(A,B,eps=1e-8){
C=A%*%B
a=trace(C)
b=eigen(C)$values
b=Re(b)
ind=which(b>eps)
b=sum(log(b[ind]))
return(a-b-ncol(A))
}

Penamatrix=function(H2){
S=diag(H2)
for(i in 1:length(H2)){
for(j in i:length(H2)){
S[i,j]=S[j,i]=1/sqrt(H2[i]*H2[j])  
}
}
diag(S)=0
return(S)
}

entropy.mcp.spearman.sampling=function(BETA,Rnoise,lamvec=c(3:12)/100,max.eps=0.001,max.iter=25,rho=0.25,mineig=0.01,subtime=100,subfrac=0.66,subthres=0.95,alpha=0.05,Const="none",variance="MAD"){
m=dim(BETA)[1];p=dim(BETA)[2]
alpha=alpha*rho
if(Const[1]=="none") Const=matrix(1,p,p)
#########################tuning parameter selection##############################  
subvecerror=matrix(0,length(lamvec),subtime)
Thetalist=array(NA,c(length(lamvec),subtime,p,p))

for(j in 1:subtime){
indsub=sample(m,round(subfrac*m),replace=F)  
BETAS=BETA[indsub,]
S=spearmancov(BETAS,method=variance)*nrow(BETAS)
M=Rnoise*dim(BETAS)[1]
S1=cov2cor1(S-M)
Theta0=glasso::glasso(S1,0.1)$wi

BETAS2=BETA[-indsub,]
S2=spearmancov(BETAS2,method=variance)*nrow(BETAS2)
M2=Rnoise*dim(BETAS2)[1]
S2=cov2cor1(S2-M2)
S2=MCPthreshold(S2,2*sqrt(log(p)/m),ga=3)
for(i in 1:length(lamvec)){
Theta=Theta0
Theta1=Theta0*0
Delta1=Theta
Gamma1=Delta1*0
Delta2=Theta
Gamma2=Delta2*0
error=norm(Theta-Theta1,"f")/sqrt(p)
iter=0
while(error>max.eps & iter<max.iter){
Theta1=Theta
Q=S1+Gamma1-rho*Delta1
Theta=(-Q+matrixsqrt(Q%*%Q+4*rho*diag(p))$w)/(2*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
if(iter %% 5==0 & min(eigen(Theta)$values)<mineig){
Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
Theta=(-Q+matrixsqrt(Q%*%Q+8*rho*diag(p))$w)/(4*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
Delta2=positveadj(Theta+Gamma2/rho,min.eps=mineig)
Gamma2=Gamma2+rho*(Theta-Delta2)
}
iter=iter+1
error=norm(Theta-Theta1,"f")/sqrt(p)
}  
df=(sum(Delta1!=0)-p)/2
subvecerror[i,j]=entropyloss(S2,Delta1)+(log(p*(p-1)/2)+log(m))/m*df  
Thetalist[i,j,,]=Delta1  
}
}

istar=which.min(rowMedian(subvecerror))
Thetalist=Thetalist[istar,,,]
K=Thetalist[1,,]*0
Thetasub=K
for(i in 1:subtime){
K=K+(Thetalist[i,,]!=0)/subtime
}

S=spearmancov(BETA,method=variance)*nrow(BETA)
M=Rnoise*dim(BETA)[1]
S1=cov2cor1(S-M)
Theta0=glasso::glasso(S1,0.1)$wi
Theta=Theta0
Theta1=Theta0*0
Delta1=Theta
Gamma1=Delta1*0
Delta2=Theta
Gamma2=Delta2*0
error=norm(Theta-Theta1,"f")/sqrt(p)
iter=0
while(error>max.eps & iter<(2*max.iter)){
Theta1=Theta
Q=S1+Gamma1-rho*Delta1
Theta=(-Q+matrixsqrt(Q%*%Q+4*rho*diag(p))$w)/(2*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const*(K>subthres)
Gamma1=Gamma1+rho*(Theta-Delta1)
if(iter %% 5==0 & min(eigen(Theta)$values)<mineig){
Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
Theta=(-Q+matrixsqrt(Q%*%Q+8*rho*diag(p))$w)/(4*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
Delta2=positveadj(Theta+Gamma2/rho,min.eps=mineig)
Gamma2=Gamma2+rho*(Theta-Delta2)
}
iter=iter+1
error=norm(Theta-Theta1,"f")/sqrt(p)
}  
Theta0=Delta1
AA=MRGNP(Theta0,Thetalist)
colnames(Theta0)=rownames(Theta0)=colnames(K)=rownames(K)=colnames(AA)=rownames(AA)=colnames(BETA)

A=list()  
A$Theta=Theta0
A$Pvalue=AA
A$K=K
A$cv.error=subvecerror
A$R=S1
return(A)
}

dtrace.mcp.spearman.sampling=function(BETA,Rnoise,lamvec=c(3:12)/100,max.eps=0.001,max.iter=25,rho=0.25,mineig=0.01,subtime=100,subfrac=0.66,subthres=0.95,alpha=0.05,Const="none",variance="MAD"){
m=dim(BETA)[1];p=dim(BETA)[2]
alpha=alpha*rho
if(Const[1]=="none") Const=matrix(1,p,p)
#########################tuning parameter selection##############################  
subvecerror=matrix(0,length(lamvec),subtime)
Thetalist=array(NA,c(length(lamvec),subtime,p,p))

for(j in 1:subtime){
  indsub=sample(m,round(subfrac*m),replace=F)  
  BETAS=BETA[indsub,]
  S=spearmancov(BETAS,method=variance)*nrow(BETAS)
  M=Rnoise*dim(BETAS)[1]
  S1=cov2cor1(S-M)
  Theta0=glasso::glasso(S1,0.1)$wi
  
  BETAS2=BETA[-indsub,]
  S2=spearmancov(BETAS2,method=variance)*nrow(BETAS2)
  M2=Rnoise*dim(BETAS2)[1]
  S2=cov2cor1(S2-M2)
  S2=MCPthreshold(S2,2*sqrt(log(p)/m),ga=3)
  for(i in 1:length(lamvec)){
    Theta=Theta0
    Theta1=Theta0*0
    Delta1=Theta
    Gamma1=Delta1*0
    Delta2=Theta
    Gamma2=Delta2*0
    error=norm(Theta-Theta1,"f")/sqrt(p)
    iter=0
    while(error>max.eps & iter<max.iter){
      Theta1=Theta
      Theta=Dtrace(S1+rho*diag(p),diag(p)+rho*Delta1-Gamma1)
      Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
      Delta1=matrix(Delta1,p,p)*Const
      Gamma1=Gamma1+rho*(Theta-Delta1)
      if(iter %% 5==0 & min(eigen(Theta)$values)<mineig){
        Theta=Dtrace(S1+2*rho*diag(p),diag(p)+rho*Delta1-Gamma1+rho*Delta2-Gamma2)
        Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
        Delta1=matrix(Delta1,p,p)*Const
        Gamma1=Gamma1+rho*(Theta-Delta1)
        Delta2=positveadj(Theta+Gamma2/rho,min.eps=mineig)
        Gamma2=Gamma2+rho*(Theta-Delta2)
      }
      iter=iter+1
      error=norm(Theta-Theta1,"f")/sqrt(p)
    }  
    df=(sum(Delta1!=0)-p)/2
    subvecerror[i,j]=entropyloss(S2,Delta1)+(log(p*(p-1)/2)+log(m))/m*df  
    Thetalist[i,j,,]=Delta1  
  }
}

istar=which.min(rowMedian(subvecerror))
Thetalist=Thetalist[istar,,,]
K=Thetalist[1,,]*0
Thetasub=K
for(i in 1:subtime){
  K=K+(Thetalist[i,,]!=0)/subtime
}

S=spearmancov(BETA,method=variance)*nrow(BETA)
M=Rnoise*dim(BETA)[1]
S1=cov2cor1(S-M)
Theta0=glasso::glasso(S1,0.1)$wi
Theta=Theta0
Theta1=Theta0*0
Delta1=Theta
Gamma1=Delta1*0
Delta2=Theta
Gamma2=Delta2*0
error=norm(Theta-Theta1,"f")/sqrt(p)
iter=0
while(error>max.eps & iter<(2*max.iter)){
  Theta1=Theta
  Theta=Dtrace(S1+rho*diag(p),diag(p)+rho*Delta1-Gamma1)
  Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
  Delta1=matrix(Delta1,p,p)*Const*(K>subthres)
  Gamma1=Gamma1+rho*(Theta-Delta1)
  if(iter %% 5==0 & min(eigen(Theta)$values)<mineig){
    Theta=Dtrace(S1+2*rho*diag(p),diag(p)+rho*Delta1-Gamma1+rho*Delta2-Gamma2)
    Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
    Delta1=matrix(Delta1,p,p)*Const
    Gamma1=Gamma1+rho*(Theta-Delta1)
    Delta2=positveadj(Theta+Gamma2/rho,min.eps=mineig)
    Gamma2=Gamma2+rho*(Theta-Delta2)
  }
  iter=iter+1
  error=norm(Theta-Theta1,"f")/sqrt(p)
}  
Theta0=Delta1
AA=MRGNP(Theta0,Thetalist)
colnames(Theta0)=rownames(Theta0)=colnames(K)=rownames(K)=colnames(AA)=rownames(AA)=colnames(BETA)

A=list()  
A$Theta=Theta0
A$Pvalue=AA
A$K=K
A$cv.error=subvecerror
A$R=S1
return(A)
}

entropy.mcp.linear.sampling=function(BETA,Rnoise,lamvec=c(3:12)/100,max.eps=0.001,max.iter=25,rho=0.25,mineig=0.01,subtime=100,subfrac=0.66,subthres=0.95,alpha=0.05,Const="none"){
m=dim(BETA)[1];p=dim(BETA)[2]
alpha=alpha*rho
if(Const[1]=="none") Const=matrix(1,p,p)
#########################tuning parameter selection##############################  
subvecerror=matrix(0,length(lamvec),subtime)
Thetalist=array(NA,c(length(lamvec),subtime,p,p))

for(j in 1:subtime){
indsub=sample(m,round(subfrac*m),replace=F)  
BETAS=BETA[indsub,]
S=t(BETAS)%*%BETAS
M=Rnoise*length(indsub)
S1=cov2cor1(S-M)
Theta0=glasso::glasso(S1,0.1)$wi
BETAS2=BETA[-indsub,]
S2=t(BETAS2)%*%BETAS2
M2=dim(BETAS2)[1]*Rnoise
S2=cov2cor1(S2-M2)
S2=MCPthreshold(S2,2*sqrt(log(p)/m),ga=3)
for(i in 1:length(lamvec)){
Theta=Theta0
Theta1=Theta0*0
Delta1=Theta
Gamma1=Delta1*0
Delta2=Theta
Gamma2=Delta2*0
error=norm(Theta-Theta1,"f")/sqrt(p)
iter=0
while(error>max.eps & iter<max.iter){
Theta1=Theta
Q=S1+Gamma1-rho*Delta1
Theta=(-Q+matrixsqrt(Q%*%Q+4*rho*diag(p))$w)/(2*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
if(iter %% 4==0 & min(eigen(Theta)$values)<mineig){
Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
Theta=(-Q+matrixsqrt(Q%*%Q+8*rho*diag(p))$w)/(4*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[i]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
Delta2=positveadj(Theta+Gamma2/rho,min.eps=mineig)
Gamma2=Gamma2+rho*(Theta-Delta2)
}
iter=iter+1
error=norm(Theta-Theta1,"f")/sqrt(p)
}  
df=(sum(Delta1!=0)+p)/2
subvecerror[i,j]=entropyloss(S2,Delta1)+(log(p*(p-1)/2)+log(m))/m*df  
Thetalist[i,j,,]=Delta1  
}
}

istar=which.min(rowMedian(subvecerror))
Thetalist=Thetalist[istar,,,]
K=Thetalist[1,,]*0
Thetasub=K
for(i in 1:subtime){
K=K+(Thetalist[i,,]!=0)/subtime
}

S=t(BETA)%*%BETA
M=Rnoise*dim(BETA)[1]
S1=cov2cor1(S-M)
Theta0=glasso::glasso(S1,0.1)$wi
Theta=Theta0
Theta1=Theta0*0
Delta1=Theta
Gamma1=Delta1*0
Delta2=Theta
Gamma2=Delta2*0
error=norm(Theta-Theta1,"f")/sqrt(p)
iter=0
while(error>max.eps & iter<(max.iter*2)){
Theta1=Theta
Q=S1+Gamma1-rho*Delta1
Theta=(-Q+matrixsqrt(Q%*%Q+4*rho*diag(p))$w)/(2*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
if(iter %% 4==0 & min(eigen(Theta)$values)<mineig){
Q=S1+Gamma1+Gamma2-rho*(Delta1+Delta2)
Theta=(-Q+matrixsqrt(Q%*%Q+8*rho*diag(p))$w)/(4*rho+alpha)
Delta1=mcp(vec(Theta+Gamma1/rho),lamvec[istar]/rho,ga=3)
Delta1=matrix(Delta1,p,p)*Const
Gamma1=Gamma1+rho*(Theta-Delta1)
Delta2=positveadj(Theta+Gamma2/rho,min.eps=mineig)
Gamma2=Gamma2+rho*(Theta-Delta2)
}
iter=iter+1
error=norm(Theta-Theta1,"f")/sqrt(p)
}  
Theta0=Delta1
AA=MRGNP(Theta0,Thetalist)
colnames(Theta0)=rownames(Theta0)=colnames(K)=rownames(K)=colnames(AA)=rownames(AA)=colnames(BETA)

A=list()  
A$Theta=Theta0
A$Pvalue=AA
A$K=K
A$cv.error=subvecerror
A$R=S1
return(A)
}

graph.susie.linear=function(BETA,Rnoise,L=10,pip.thres=0.5){
  n=dim(BETA)[1];p=dim(BETA)[2]
  S=t(BETA)%*%BETA
  M=Rnoise*n
  S1=cov2cor1(S-M)
  Theta=S1*0
  for(i in 1:p){
    W11=S1[-i,-i]
    s1=c(S1[i,-i])
    fit=susie_suff_stat(XtX=W11*n,Xty=s1*n,yty=n*S1[i,i],n=n,L=L)
    beta=coef(fit)[-1]*(fit$pip>pip.thres)
    w12=c(W11%*%beta)
    theta22=1/(S1[i,i]-sum(w12*beta))
    theta12=-beta*theta22
    Theta[i,i]=theta22
    Theta[-i,i]=Theta[i,-i]=theta12
  }
  Theta=as.symmeter(Theta)
  rownames(Theta)=colnames(Theta)=colnames(BETA)
  return(A=list(Theta=Theta))
}

graph.susie.spearman=function(BETA,Rnoise,L=10,pip.thres=0.5){
  n=dim(BETA)[1];p=dim(BETA)[2]
  S=spearmancov(BETA,method="MAD")*nrow(BETA)
  M=Rnoise*n
  S1=cov2cor1(S-M)
  Theta=S1*0
  for(i in 1:p){
      W11=S1[-i,-i]
      s1=c(S1[i,-i])
      fit=susie_suff_stat(XtX=W11*n,Xty=s1*n,yty=n*S1[i,i],n=n,L=L)
      beta=coef(fit)[-1]*(fit$pip>pip.thres)
      w12=c(W11%*%beta)
      theta22=1/(S1[i,i]-sum(w12*beta))
      theta12=-beta*theta22
      Theta[i,i]=theta22
      Theta[-i,i]=Theta[i,-i]=theta12
  }
  Theta=as.symmeter(Theta)
  rownames(Theta)=colnames(Theta)=colnames(BETA)
  return(A=list(Theta=Theta))
}

graph.susie.kendall=function(BETA,Rnoise,L=10,pip.thres=0.5){
  n=dim(BETA)[1];p=dim(BETA)[2]
  S=kendallcov(BETA,method="MAD")*nrow(BETA)
  M=Rnoise*n
  S1=cov2cor1(S-M)
  Theta=S1*0
  for(i in 1:p){
    W11=S1[-i,-i]
    s1=c(S1[i,-i])
    fit=susie_suff_stat(XtX=W11*n,Xty=s1*n,yty=n*S1[i,i],n=n,L=L)
    beta=coef(fit)[-1]*(fit$pip>pip.thres)
    w12=c(W11%*%beta)
    theta22=1/(S1[i,i]-sum(w12*beta))
    theta12=-beta*theta22
    Theta[i,i]=theta22
    Theta[-i,i]=Theta[i,-i]=theta12
  }
  Theta=as.symmeter(Theta)
  rownames(Theta)=colnames(Theta)=colnames(BETA)
  return(A=list(Theta=Theta))
}

