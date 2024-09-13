vec=function(a){
  as.vector(a)
}

MRGNP=function(A,B){
p=dim(A)[2];P=diag(p)
for(i in 1:p){
  for(j in i:p){
  a=A[i,j];b=sd(B[,i,j])
  P[i,j]=P[j,i]=pchisq((a/b)^2,1,lower.tail=F)
  }
}
return(P)
}

iqr=function(x){
s=(quantile(x,0.75)-quantile(x,0.25))/1.349
return(s)
}

apd=function(x){
p=length(x)
D=diag(p) 
for(i in 1:p){
D[,i]=abs(x-x[i])
}
d=colMedian(D)
s=median(d)*1.1926
return(s)
}

as.symmeter <- function(matrix) {
  n <- nrow(matrix)
  if (n != ncol(matrix)) {
    stop("The matrix must be square.")
  }
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (abs(matrix[i, j]) < abs(matrix[j, i])) {
        matrix[j, i] <- matrix[i, j]
      } else {
        matrix[i, j] <- matrix[j, i]
      }
    }
  }
  
  return(matrix)
}

normtransform=function(X){
  X1=apply(X,2,rank)/length(X[,1])
  X1=qnorm(X1)
  X1[abs(X1)>5]=5*sign(X1[abs(X1)>5])
  return(X1)
}

climeselect=function(X,climefit){
p=length(climefit$lambda)
bicvec=c(1:p)
R=cov2cor(matrixMultiply(t(X),X))
for(i in 1:p){
A=climefit$Omegalist[[i]]
A[abs(A)<1e-4]=0
climefit$Omegalist[[i]]=A
df=(sum(A!=0)+ncol(A))/2
bicvec[i]=entropyloss(R,A)+log(nrow(X))/nrow(X)*df
}
climefit$bic=bicvec
climefit$Theta=climefit$Omegalist[[which.min(bicvec)]]
return(climefit)
}

thresholding=function(A,lambda){
  A1=cov2cor(A)
  d=sqrt(diag(A))
  A2=abs(A1)
  A2=A2-lambda
  A2[A2<0]=0
  diag(A2)=1
  A2=A2*sign(A1)
  A=t(t(A2*d)*d)
  return(A)
}

colSDmed=function(A){
  B=A
  p=ncol(A)
  for(i in 1:p){
    ind=which(A[,i]!=0)
    B[,i]=sd(A[ind,i])
  }
  return(B)
}

biggwas=function(x,G){
  x=as.vector(x)
  ux=mean(x)
  vx=var(x);vx=as.numeric(vx)
  ug=colMeans(G)
  G=t(t(G)-ug)
  vg=colSums(G^2)
  b=(t(G)%*%(x-ux))/vg
  sdb=(vx-b^2*vg/length(x))/length(x)
  A=list()
  A$est=as.vector(b)
  A$std=as.vector(sqrt(sdb))
  return(A)
}

mbiggwas=function(X,G1,N){
  p=dim(X)[2]
  m=dim(G1)[2]
  B=matrix(0,m,p)
  SDB=B
  for(i in 1:p){
    x=X[N[,i],i]
    G=G1[N[,i],]
    x=as.vector(x)
    ux=mean(x)
    vx=var(x);vx=as.numeric(vx)
    ug=colMeans(G)
    G=t(t(G)-ug)
    vg=colSums(G^2)
    b=(t(G)%*%(x-ux))/vg
    b=as.vector(b)
    sdb=(vx-b^2*vg/length(x))/length(x)
    B[,i]=b
    SDB[,i]=sqrt(as.vector(sdb))
  }
  A=list()
  A$est=B
  A$std=SDB
  return(A)
}

colMedian=function(A){
a=apply(A,2,median)
return(a)
}

rowMedian=function(A){
  a=apply(A,1,median)
  return(a)
}

vary=function(Vg,theta,rho,hy){
  b=2*rho*theta  
  c=theta^2-Vg*(1/hy-1)
  d=-b+sqrt(b^2-4*c)
  d=d/2
  d=as.numeric(d)
  return(d^2)
}

mvary=function(Vg,theta,SigmaUU,SigmaUV,hy){
  b=2*sum(theta*SigmaUV)
  c=Vg*(1-1/hy)+t(theta)%*%SigmaUU%*%theta
  d=-b+sqrt(max(b^2-4*c,0))
  d=d/2
  d=as.numeric(d)
  d=max(d,0)
  return(d^2)
}

svar=function(x){
  x=as.vector(x)
  mx=mean(x)
  v=sum((x-mx)^2)/length(x[-1])  
}

blockdiag=function(A,B){
  if(length(A)>1){
    p1=dim(A)[2];p2=dim(B)[2]
    C=diag(p1+p2)
    C[1:p1,1:p1]=A
    C[c((p1+1):(p1+p2)),c((p1+1):(p1+p2))]=B
  }else{
    p1=1; 
    p2=dim(B)[2]
    C=diag(p1+p2)
    C[1,1]=A
    C[c(2:(1+p2)),c(2:(1+p2))]=B
  }
  return(C)
}

CScov=function(p,rho){
  A=diag(p)*(1-rho)+matrix(rho,p,p)
  return(A)
}

ARcov=function(p,rho){
  s=c(1:p)
  for(i in 1:p) s[i]=rho^(i-1)
  return(toeplitz(s))
}

innercov=function(hata,hatb,theta,byse,bXse,rxy){
  w=cbind(1,hatb)
  err=hata-w%*%theta
  m=dim(w)[1];p=dim(w)[2]
  G=matrix(0,m,p)
  for(i in 1:m){
    r=err[i]*w[i,]-bXse[i]^2*theta+rxy*byse[i]*bXse[i]
    G[i,]=r  
  }
  return(G)
}

.usumSuu=function(bXse,Rxx){
  n=length(bXse)
  Suu=0
  for(i in 1:n){
    Suui=bXse[i]^2*Rxx
    Suu=Suu+Suui
  }
  return(Suu)
}

spearmancov=function(A,method="MAD"){
  p=dim(A)[2]
  s=c(1:p)
  if(method=="MAD"){
  for(i in 1:p){
    s[i]=median(abs(median(A[,i])-A[,i]))*1.483
  }}
  if(method=="MedianOfMean"){
  s=sqrt(medianofmean(A))
  }
  R=cor(A,method="spearman")
  R=2*sin(R*pi/6)
  S=diag(s)%*%R%*%diag(s)
  return(S)
}

kendallcov=function(A,method="MAD"){
  p=dim(A)[2]
  s=c(1:p)
  if(method=="MAD"){
    for(i in 1:p){
      s[i]=median(abs(median(A[,i])-A[,i]))*1.483
    }}
  if(method=="MedianOfMean"){
    s=sqrt(medianofmean(A))
  }
  R=cor(A,method="kendall")
  R=sin(R*pi/2)
  S=diag(s)%*%R%*%diag(s)
  return(S)
}

.usumSuv=function(bXse,byse,rxy){
  n=length(bXse)
  Suv=0
  for(i in 1:n){
    Suvi=bXse[i]*rxy*byse[i]
    Suv=Suv+Suvi
  }
  return(Suv)
}


rowdotproduct=function(A,a){
  a=as.vector(a)
  p=length(A[1,])
  for(i in 1:p){
    A[,i]=A[,i]*a
  }
  return(A)
}

cov2cor1=function(A,kappa=100){
fit=matrixEigen(A)
d=c(fit$values)
eps=max(d)/kappa
d[d<eps]=eps
B=matrixMultiply(fit$vectors,t(fit$vectors)*d)
colnames(B)=rownames(B)=colnames(A)
return(cov2cor(B))
}

.sumSuu=function(bXse,Rxx){
  n=length(bXse[,1])
  Suu=0
  for(i in 1:n){
    Suui=diag(bXse[i,])%*%Rxx%*%diag(bXse[i,])
    Suu=Suu+Suui
  }
  return(Suu)
}

.sumSuv=function(bXse,byse,rxy){
  n=length(bXse[,1])
  Suv=0
  for(i in 1:n){
    Suvi=bXse[i,]*rxy*byse[i]
    Suv=Suv+Suvi
  }
  return(Suv)
}

mvvarv=function(Vbb,Vuu,Ruv,theta,hy){
  Vu=diag(Vuu)
  b=2*sum(theta*Ruv*sqrt(Vu))
  c=t(theta)%*%(Vbb+Vuu)%*%theta-t(theta)%*%Vbb%*%theta/hy
  vv=(-b+sqrt(b^2-4*c))/2
  return(vv^2)
}


dSCAD=function(a,lam,gamma=3.7){
  a=abs(a)
  z=a
  z[a<lam]=lam
  z[a>lam]=(gamma*lam-z[a>lam])/(gamma-1)
  z[a>(gamma*lam)]=0  
  return(z)
}

dMT=function(a,lam,sig=3.7){
  r=abs(a/lam)
  z=exp(-r/sig)*(1+r/sig)
  z=z*lam
  return(z)
}


dMCP=function(a,lam,gamma=3.7){
  a=abs(a)
  z=lam-a/gamma
  z[a>(gamma*lam)]=0  
  return(z)
}

soft=function(a,b){
  c=abs(a)-b
  c[c<0]=0
  c=c*sign(a)
  return(c)
}

positveadj=function(A,min.eps=0.001){
a=matrixEigen(A)
d=c(a$values)
d[d<min.eps]=0
B=matrixMultiply(a$vectors,t(a$vectors)*d)
return(B)
}


standard=function(b){
  p=ncol(b)
  for(j in 1:p){
    s=b[,j]
    s=(s-mean(s))/sd(s)
    b[,j]=s
  }
  return(b)
}

Uerrorvar=function(Sbb,hxvec){
  p=length(hxvec)
  Vb=d=diag(Sbb)
  D=Vb/hxvec
  return(D)
}

verrorvar=function(Sbb,Ruv,Du,theta,hy){
  p=length(theta)
  Suu=diag(sqrt(Du))%*%Ruv[1:p,1:p]%*%diag(sqrt(Du));
  ruv=Ruv[1:p,p+1]
  svv=mvvarv(Vbb=Sbb,Vuu=Suu,Ruv=ruv,theta=theta,hy=hy)  
  return(svv)
}

CovXy=function(Sbb,Suv,Rover,theta){
  p=length(theta)
  Suu=Suv[1:p,1:p]
  Svv=Suv[p+1,p+1]
  Suv=Suv[1:p,p+1]
  Sxy=diag(p+1)
  Sxy[1:p,1:p]=Suu+Sbb
  Sxy[1:p,p+1]=Sxy[p+1,1:p]=(Suu+Sbb)%*%theta+Suv
  Sxy[1+p,1+p]=t(theta)%*%(Suu+Sbb)%*%theta+2*sum(theta*Suv)+Svv
  Vxy=Sxy*Rover 
  return(Vxy)
}

ARmatrix=function(rho,p){
  l=c(1:p)
  for(i in 2:p){
    l[i]=rho^(i-1)
  }
  L=toeplitz(l)
  return(L)
}

CSmatrix=function(rho,p){
  L=diag(p)*(1-rho)+matrix(rho,p,p)
  return(L)
}

EXPmatrix=function(rho,d){
D=diag(d)
p=length(d)
for(i in 1:p){
  for(j in i:p){
    s=abs(d[i]-d[j])
    D[i,j]=D[j,i]=exp(-s/rho)
  }
}
return(D)
}

SPHmatrix=function(rho,d){
p=length(d)
D=diag(p)
for(j in 1:p){
s=abs(d-d[j])
g=1-1.5*s/rho+0.5*(s/rho)^3
g=as.vector(g)
g[s>rho]=0
D[,j]=g
}
return(D)
}

corrmatrix=function(A){
  a=1/sqrt(diag(A))
  B=diag(a)%*%A%*%diag(a)
  return(B)
}

SCADthreshold=function(S,lam,k=3){
  d=sqrt(diag(S))
  R=cov2cor(S)
  s=as.vector(R)
  s1=abs(s)
  s2=scad(s1,lam)
  s2=s2*sign(s)
  S1=matrix(s2,ncol(S),ncol(S))
  diag(S1)=1
  S2=t(S1*d)*d
  return(S2)
}

MCPthreshold=function(S,lam,ga=3){
  d=sqrt(diag(S))
  R=cov2cor(S)
  s=as.vector(R)
  s1=abs(s)
  s2=mcp(s1,lam,ga)
  s2=s2*sign(s)
  S1=matrix(s2,ncol(S),ncol(S))
  diag(S1)=1
  S2=diag(d)%*%S1%*%diag(d)
  return(S2)
}

tapering=function(A,k){
  p=ncol(A)
  l=rep(1,p)*0
  l[1:k]=1
  for(i in (k+1):min((2*k),p)){
    l[i]=2-i/k
  }
  return(toeplitz(l)*A)
}


scad=function(a,lam,ga=3.7){
b=abs(a)
z=soft(a,lam)
z[which(b>(2*lam)&b<=(ga*lam))]=soft(a[which(b>(2*lam)&b<=(ga*lam))],ga*lam/(ga-1))/(1-1/(ga-1))
z[which(b>(ga*lam))]=a[which(b>(ga*lam))]
return(z)
}

mcp=function(a,lam,ga=3.7){
b=abs(a)
z=soft(a,lam)/(1-1/ga)
z[which(b>(ga*lam))]=a[which(b>(ga*lam))]
return(z)
}

var.error=function(byse,bXse,theta,Rxx,rxy){
  m=length(by)
  vr=c(1:m)*0
  for(ss in 1:m){
    Suui=diag(bXse[ss,])%*%Rxx%*%diag(bXse[ss,])
    Suvi=byse[ss]*bXse[ss,]*rxy
    vr[ss]=byse[ss]^2+t(theta)%*%Suui%*%theta-2*sum(theta*Suvi)
  }
  return(vr)
}

columsd=function(A){
n=dim(A)[2]
s=c(1:n)
for(i in 1:n){
s[i]=sd(A[,i])
}
return(s)
}

mad=function(a){
b=1.483*median(abs(a-median(a)))
return(b)
}

inverseARcov=function(p,rho){
a=c(1+rho^2,-rho,rep(0,p-2)) 
A=toeplitz(a)  
A[1,1]=A[p,p]=1  
A=A/(1-rho^2)
A=Matrix::Matrix(A,sparse=T)
return(A)
}

inverseCScov=function(p,rho){
k1=(p-1)*rho^2-(p-2)*rho-1
a0=-((p-2)*rho+1)/k1
a1=rho/k1
A=a0*diag(p)+a1*(matrix(1,p,p)-diag(p))
return(A)
}

trace=function(A){
a=sum(diag(A))
return(a)
}

normvec=function(a){
return(sqrt(sum(a^2)))
}

coverage=function(a,b,method="BC"){
ind=which(a!=0)
pv=1+a*0
thres1=0.05/length(ind)
pv[ind]=1-pchisq(a[ind]^2/b[ind]^2,1)
issign1=as.numeric(pv<thres1)
if(method=="BC"){
return(issign1)
}
if(method!="BC"){
  pv1=p.fdr(pvalues=pv)
  issign2=as.numeric(pv1$fdrs<0.05,adjust.method=method)
  return(issign2)
}
}

tptn=function(thetahat,theta){
b=as.numeric(theta!=0) 
a=as.numeric(thetahat!=0)
ind=which(theta!=0)
tp=sum(a[ind]!=b[ind])==0
tn=sum(a[-ind]!=b[-ind])==0
A=list(tp=as.numeric(tp),tn=as.numeric(tn))
return(A)
}

infer=function(by,bX,theta){
ind=which(theta!=0)
r=by-bX%*%theta
see=mad(r)
B=diag(ncol(bX))*0
if(sum(ind)>0){
B[ind,ind]=see*solve(t(bX[,ind])%*%bX[,ind])
}
return(B)
}

Zmatrix=function(SNP,GeneSymbol,Zscore){
snp=unique(SNP)
gene=unique(GeneSymbol)
m=length(snp);p=length(gene)
Z=matrix(NA,m,p)
for(i in 1:m){
  for(j in 1:p){
    ind=which(SNP==snp[i]&GeneSymbol==gene[j])
    if(length(ind)==1){Z[i,j]=Zscore[ind]}
    if(length(ind)>1){Z[i,j]=Zscore[ind[1]]}
  }
}
colnames(Z)=gene
rownames(Z)=snp
return(Z)
}

demissing=function(M,rowthres=0.2,colthres=0.2,rowfirst=T){
M1=is.na(M)
if(rowfirst==T){
M2=rowMeans(M1)
ind1=which(M2<=rowthres)
M2=colMeans(M1[ind1,])
ind2=which(M2<=colthres)
}else{
M2=colMeans(M1)
ind2=which(M2<=colthres) 
M2=rowMeans(M1[,ind2])
ind1=which(M2<=rowthres)
}  
M3=M[ind1,ind2]
return(M3)
}


allele_harmonise <- function(ref_panel, gwas_data) {
  # Make sure the reference panel has columns SNP, A1 and A2
  stopifnot(all(c("SNP", "A1", "A2") %in% colnames(ref_panel)))
  # Make sure the GWAS data has columns SNP, A1, A2, and Tstat
  stopifnot(all(c("SNP", "A1", "A2", "Zscore") %in% colnames(gwas_data)))
  
  # Merge the reference panel with the GWAS data
  merged_data <- merge(gwas_data, ref_panel, by = "SNP")
  
  # Remove any SNPs where A1 and A2 do not match between the reference panel and the GWAS data
  ind=which((toupper(merged_data$A1.x) == toupper(merged_data$A1.y) & toupper(merged_data$A2.x) == toupper(merged_data$A2.y))
            |(toupper(merged_data$A1.x) == toupper(merged_data$A2.y) & toupper(merged_data$A2.x) == toupper(merged_data$A1.y)))
  merged_data=merged_data[ind,]
  # Multiply Tstat by -1 if A1 and A2 are opposed between the reference panel and the GWAS data
  merged_data$Zscore <- ifelse(toupper(merged_data$A1.x) == toupper(merged_data$A2.y) & toupper(merged_data$A2.x) == toupper(merged_data$A1.y), -1 * merged_data$Zscore, merged_data$Zscore)
  
  colnames(merged_data)[which(names(merged_data) == "A1.y")] <- "A1"
  colnames(merged_data)[which(names(merged_data) == "A2.y")] <- "A2"
  colnames(merged_data)[which(names(merged_data) == "A1.x")] <- "A11"
  colnames(merged_data)[which(names(merged_data) == "A2.x")] <- "A21"
  # Convert A1 and A2 to uppercase
  merged_data$A1 <- toupper(merged_data$A1)
  merged_data$A2 <- toupper(merged_data$A2)
  
  return(merged_data)
}

bimin=function(mat){
    min_element <- min(mat)
    min_indices <- which(mat == min_element, arr.ind = TRUE)
    if (nrow(min_indices) > 1) {
      min_indices <- min_indices[nrow(min_indices), ]
    }
    return(min_indices)
}

lowrankadj=function(A,eps){
fit=matrixEigen(A)
d=c(fit$value)
d1=cumsum(d)/sum(d)
ind=min(which(d1>eps))
if(ind<ncol(A)){
d[ind:ncol(A)]=0  
}
B=matrixMultiply(fit$vector,t(fit$vector)*d)
return(B)
}

matrixsqrt=function(A){
fit=matrixEigen(A)
d=c(fit$value)  
d[d<0]=0
d1=d*0
d1[d>0]=1/d[d>0]
d=sqrt(d)
d1=sqrt(d1)
A=matrixMultiply(fit$vector,t(fit$vector)*d)
B=matrixMultiply(fit$vector,t(fit$vector)*d1)
C=list(w=A,wi=B)
return(C)
}

distancematrix=function(b){
p=length(b)
a=matrix(0,length(b),length(b))  
for(i in 1:p){
a[,i]=abs(b[i]-b)
}
return(a)  
}

pratt.index=function(y,x,w,z){
  y=as.vector(y);x=as.vector(x);w=as.vector(w);z=as.vector(z)
  B=cbind(x,w,z)  
  s=c(sum(abs(x))>0,sum(abs(w)>0),sum(abs(z)>0))
  if(sum(s==3)){
    fit=lm(y~x+w+z-1)  
    theta=fit$coefficient  
    delta=cor(y,B)
    para=theta*delta
  }
  if(sum(s<2)&sum(s)>0){
    ind=which(s>0)
    fit=lm(y~B[,ind]-1) 
    theta=delta=B[1,]*0
    theta[ind]=fit$coefficient
    varb=colMeans(B^2)
    delta[ind]=cov(y,B[,ind])/varb[ind]
    para=theta*delta
  }
  if(sum(s)==0) para=0*B[1,]  
  return(para)  
}

diagcomb=function(A,B){
p1=ncol(A);p2=ncol(B)
C=diag(p1+p2)
C[1:p1,1:p1]=A
C[-c(1:p1),-c(1:p1)]=B
return(C)
}

lowrankimpute=function(X,R,K=round(ncol(X)/2),iter=10){
X1=X;indna=is.na(X)
Theta=solve(R)
X1[indna]=0
for(i in 1:iter){
eig=matrixEigen(t(X1)%*%Theta%*%X1)
Gamma=eig$vector[,1:K]
L=X1%*%Gamma%*%t(Gamma)
X1[indna]=L[indna]
}
return(X1)
}

symm=function(A){
p=ncol(A)
B=A
for(i in 1:p){
for(j in i:p){
a=A[i,j];b=A[j,i]
s=which.min(abs(c(a,b)))
if(s==1){
B[i,j]=B[j,i]=a
}else{
B[i,j]=B[j,i]=b  
}
}
}
return(B)  
}

directed=function(A){
  p=ncol(A)
  B=A
  for(i in 1:p){
    for(j in i:p){
      a=A[i,j];b=A[j,i]
      s=which.max(abs(c(a,b)))
      if(s==1){
      B[j,i]=0
      }else{
      B[i,j]=0
      }
    }
  }
  return(B)  
}

MAcov=function(p,rhovec){
  s=c(1:p)*0
  s[1:length(rhovec)]=rhovec
  return(toeplitz(s))
}

CScov=function(p,rho){
  A=matrix(rho,p,p)+(1-rho)*diag(p)
  return(A)
}

trace=function(A){
  a=sum(diag(A))
  return(a)
}

TPTN=function(A,B){
  p=ncol(A)
  P1=A!=0
  P1=matrix(as.numeric(P1),p,p)
  N1=A==0
  N1=matrix(as.numeric(N1),p,p)
  P2=B!=0
  P2=matrix(as.numeric(P2),p,p)
  N2=B==0
  N2=matrix(as.numeric(N2),p,p)  
  P2=P1*P2  
  N2=N2*N1
  tp=sum(P2-P1)
  tn=sum(N2-N1)
  return(c(tp,tn))
}

dtraceloss=function(A,B){
  a=0.5*trace(A%*%A%*%B)-trace(A)
  return(a)
}

quadraticloss=function(A,B){
  C=A%*%B-diag(ncol(A))
  a=trace(C%*%C)
  return(a)
}

LiNGAM=function(Theta,kk=F){
p=ncol(Theta)
D=diag(Theta)
B=diag(p)-diag(1/D)%*%Theta
diag(B)=0
C=abs(B)
nam=colnames(Theta)

if(kk!=F){
prep=nam[kk]
indprep=which(colnames(C)%in%prep)
C[indprep,]=0
C1=C[-indprep,-indprep]
while(length(prep)<(p-1)){
s=rowMeans(C1)
ind=which(s==0)
while(sum(ind)==0){
C1[C1==0]=100000
xx=which.min(C1)
C1[C1==100000]=0
C1[xx]=0
s=rowMeans(C1)
ind=which(s==0)
}
C[-indprep,-indprep]=C1
prep=c(colnames(C1)[ind],prep)
indprep=which(colnames(C)%in%prep)
C1=C[-indprep,-indprep]
}
prep=c(setdiff(nam,prep),prep)
}
if(kk==F){
  prep=nam[which.min(D)]
  indprep=which(colnames(C)%in%prep)
  C[indprep,]=0
  C1=C[-indprep,-indprep]
  while(length(prep)<(p-1)){
    s=rowMeans(C1)
    ind=which(s==0)
    while(sum(ind)==0){
      C1[C1==0]=100000
      xx=which.min(C1)
      C1[C1==100000]=0
      C1[xx]=0
      s=rowMeans(C1)
      ind=which(s==0)
    }
    C[-indprep,-indprep]=C1
    prep=c(colnames(C1)[ind],prep)
    indprep=which(colnames(C)%in%prep)
    C1=C[-indprep,-indprep]
  }
  prep=c(setdiff(nam,prep),prep)
}
B=B*(C!=0)
rownames(B)=colnames(B)
B=B[nam,nam]
return(B)
}