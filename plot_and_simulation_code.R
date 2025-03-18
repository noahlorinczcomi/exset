rm(list=ls(all=TRUE))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# distribution of angle of random vector in 2d space ####
library(plotly)
ar1=function(n,rho) rho^toeplitz(0:(n-1))
tr=function(x) sum(diag(x))
rho=0.9
R=ar1(2,rho)
lam=eigen(R)$values
x=seq(1e-1,2,0.1)
y=seq(1e-1,2,0.1)
fz=function(x,y) 1/(2*pi)/sqrt(prod(lam)*x*y)*exp(-1/2*(x/lam[1]+y/lam[2]))
# pracma::quad2d(fz,xa=0,ya=0,xb=100,yb=100,n=1000)
Z=outer(x,y,fz)
plot_ly(x=~x,y=~y,z=~Z,
        colors=c("#78B7C5", "#EBCC2A", "#F21A00"),
        contours=list(z=list(show=TRUE,color="black",width=2))) %>% 
  add_surface()
# randomly generate from the joint weighted independent chi-square distribution
rz=function(n,lam,ylim=10,dlim=fz(0.1,0.1)) {
  # random sampling based in geometry
  data=matrix(nr=n,nc=2)
  k=0
  while(k<n) {
    y1=runif(1,0,ylim)
    y2=runif(1,0,ylim)
    r2=runif(1,0,dlim)
    if(fz(y1,y2)>r2) {
      k=k+1
      data[k,]=c(y1,y2)
    }
  }
  data
}
sample=rz(10000,lam)
# find angles of each sampled vector
theta=apply(sample,1,function(.) atan(.[1]/.[2]))
hist(theta,col='gray90',border='gray90',pr=T,breaks=30,xlab=expression(xi))
sec2=function(x) 1/cos(x)^2
ftheta=function(x,lambdax,lambday) {
  d=lambdax/lambday
  sec2(x)/(d*pi)/(sqrt(tan(x)/d))/(1+tan(x)/d)
}
curve(ftheta(x,lam[1],lam[2]),0,pi/2,col='blue',add=T)
integrate(ftheta,0,pi/2,lam[1],lam[2])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# simulation to show controlled Type I error ####
ar1=function(n,rho) rho^toeplitz(0:(n-1))
tr=function(x) sum(diag(x))
simdata=function(niter=100000,R) {
  m=nrow(R)
  Z=mvnfast::rmvn(niter,rep(0,m),R)
  S=rowSums(Z^2)
  S
}
mc_integrate=function(lam,niter=10000,alpha=0.05) {
  p=length(lam)
  X=matrix(rchisq(niter*p,1),niter,p)
  X=X*matrix(lam,niter,p,byrow=TRUE)
  r=rowSums(X) # half plane
  tau=quantile(r,prob=1-alpha)
  tau
}
ms=c(2,50,100)
# rhos=c(-0.9,-0.5,0,0.5,0.9)
rhos=seq(-0.9,0.9,length.out=20)
alpha=0.05
niter=1e6
TAU0=TAU=RES=matrix(nr=length(ms),nc=length(rhos))
for(i in 1:length(ms)) {
  for(j in 1:length(rhos)) {
    R=ar1(ms[i],rhos[j])
    lam=eigen(R)$values
    S=simdata(niter,R)
    tau0=quantile(S,prob=1-alpha)
    tau=mc_integrate(lam,niter,alpha)
    TAU0[i,j]=tau0
    TAU[i,j]=tau
    RES[i,j]=mean(S>tau)
  }
}
pdf=data.frame(
  m=rep(ms,length(rhos)),
  rho=rep(rhos,each=length(ms)),
  type1=c(RES),
  tau0=c(TAU0),
  tau=c(TAU)
)
library(ggplot2)
p1=ggplot(pdf,aes(x=rho,y=type1,fill=factor(m))) +
  geom_bar(stat='identity',color=NA,position=position_dodge(1/length(rhos)),width=1/20) +
  geom_hline(yintercept=0.05) +
  scale_fill_manual('number of SNPs tested',values=c("#78B7C5", "#EBCC2A", "#F21A00")) +
  theme_bw() +
  lims(y=c(0,0.1)) +
  theme(legend.position='bottom',
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  labs(x='correlation between neighboring SNPs (AR1 structure)',
       y='achieved Type 1 error (0.05 target)')
p2=ggplot(pdf,aes(rho,tau,fill=factor(m))) +
  geom_abline(intercept=0,slope=1,col='gray80') +
  geom_line() +
  geom_point(pch=21,size=3) +
  scale_fill_manual('number of SNPs tested',values=c("#78B7C5", "#EBCC2A", "#F21A00")) +
  theme_bw() +
  theme(legend.position='bottom',
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank()) +
  labs(x='correlation between neighboring SNPs (AR1 structure)',
       y=expression(tau*' (95% critical value)'))
ggpubr::ggarrange(p1,p2,common.legend=T,nrow=1,ncol=2)