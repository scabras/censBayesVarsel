ldinvgamma=function (x, shape, scale = 1) 
{
  if (shape <= 0 | scale <= 0) {
    stop("Shape or scale parameter negative in dinvgamma().\n")
  }
  alpha <- shape
  beta <- scale
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
                                                        1) * log(x) - (beta/x)
  return(log.density)
}

# Log-Normal posterior with prior pi(sigma) \propto 1/sigma: (parametrizado en log(\sigma)):
log.post.lnormal.nc=function (parms, Data,b,prior=1,XtX.inv=NULL,nunc=NULL){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f= -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

log.post.lnormal=cmpfun(log.post.lnormal.nc)


##Function to calculate m1 for a set of covariates:
##Approximation marginal 1 under an alternative model, for g=1:
marg1.allpara.gfixed=function(y,rel,cens,X,k.star,nsim.marg=1000,nsim.post=5000){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=parms.hat1,fn=log.post.lnormal,control=list(fnscale=-1),Data=data1.aux,method=optimmethod,hessian=TRUE)
  var.prop=-solve(tt1$hessian)
  
  parms.hat1=tt1$par
  
  proposal=list(var=var.prop,scale=1)
  #For g=1
  res.1=simul.post.rw.gfixed(log.post.lnormal.our.gfixed,proposal,start=parms.hat1,g=1,m=nsimul,Data=as.matrix(data1.aux),cens=cens)
  res.1$par=res.1$par[-(1:(nsimul*0.3)),]
  
  mu.prop=apply(res.1$par,2,mean)
  Sigma.prop=cov(res.1$par)
  
  #With g=1
  thetas=rmvt(nsim.marg, delta = mu.prop, sigma =Sigma.prop,df=3)
  lnum=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x) log.post.lnormal.our.gfixed(thetas[x,],g=1,Data=as.matrix(data1.aux),cens=cens))
  lnum=lnum+k.star 
  lden=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x){ dmvt(thetas[x,],delta=mu.prop, sigma =Sigma.prop,df=3,log=TRUE)})
  
  marg1.comp=mean(exp(lnum-lden),na.rm=TRUE)
  return(list(marg=marg1.comp,par.post=res.1$par))
}

#Same than before but using the inverse-gamma as prior distribution over g:
marg1.allpara.igamma=function(y,rel,cens,X,k.star,nsim.marg=1000,nsim.post=5000){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=parms.hat1,fn=log.post.lnormal,control=list(fnscale=-1),Data=data1.aux,method=optimmethod,hessian=TRUE)
  var.prop=-solve(tt1$hessian)
  
  parms.hat1=tt1$par
  
  proposal=list(var=var.prop,scale=0.5)
 
   #Para la inverse-gamma:
  res.1=simul.post.rw(log.post.lnormal.our,proposal,start=parms.hat1,g.start=1,nsimul,Data=as.matrix(data1.aux),cens=cens)
  
  res.1$par=res.1$par[-(1:(nsimul*0.3)),]
  
  n.par=dim(res.1$par)[2]
  mu.prop=apply(res.1$par[,-n.par],2,mean)
  Sigma.prop=cov(res.1$par[,-n.par])
  
  #Usando g sim Inv-gamma(0.5,0.5)
  thetas=rmvt(nsim.marg, delta = mu.prop, sigma =Sigma.prop,df=3)
  gs=rinvgamma(nsim.marg,0.5,0.5)
  
  
  lnum=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x) log.post.lnormal.our(thetas[x,],g=gs[x],Data=as.matrix(data1.aux),cens=cens))
  lnum=lnum+k.star 
  lden=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x){ dmvt(thetas[x,],delta=mu.prop, sigma =Sigma.prop,df=3,log=TRUE) + ldinvgamma(gs[x],0.5,0.5)})
  
  
  marg1.comp=mean(exp(lnum-lden),na.rm=TRUE)
  return(list(marg=marg1.comp,par.post=res.1$par))
}

#Using g fixed, g=1.
#Considering the new var-covar matrix, taken into account the "quasi orthogonalization". 06/04/2017
log.prior.our.gfixed=function(parms,g=1,cens,X){
  np=length(parms)
  res=0
  n=length(cens)
  zi0=(cens-parms[2])/exp(parms[1])
  pzi0=pnorm(zi0)
  Delta0=diag(pzi0,nrow=n)
  hzi0=exp(dnorm(zi0,log=TRUE)-pnorm(zi0,lower.tail=FALSE,log.p=TRUE))
  wi=pzi0+dnorm(zi0)*(hzi0-zi0)
  ne=sum(wi)
  W=diag(wi)
  Uno=matrix(1,nrow=n,ncol=1)
  Id=diag(1,n)
  ne.inv=Uno%*%t(Uno)/ne
  W.root=diag(sqrt(wi))
  Weig=W.root%*%(Id-W.root%*%ne.inv%*%W.root)%*%W.root
  Xt.Weig.X=t(X)%*%Weig%*%X
  
  if(det(Xt.Weig.X)<1e-10 || is.na(det(Xt.Weig.X))){
    X.sca=t(X)%*%(Id-Uno%*%t(Uno)/n)%*%X
    Xt.Weig.X.inv=n*solve(X.sca)
  }else{
    Xt.Weig.X.inv=ne*solve(Xt.Weig.X)
  }
  g=1
  
  if(np>2) res=dmnorm(x=parms[-c(1,2)],mean=0,
                      varcov=exp(2*parms[1])*Xt.Weig.X.inv,log=TRUE)
  return(res)
}


#It returns the prior distribution density, for g=1 or using the inverse-gamma (comment or uncomment in any case).
log.prior.our=function(parms,g,cens,X){
  np=length(parms)
  res=0
  n=length(cens)
  zi0=(cens-parms[2])/exp(parms[1])
  pzi0=pnorm(zi0)
  Delta0=diag(pzi0,nrow=n)
  hzi0=exp(dnorm(zi0,log=TRUE)-pnorm(zi0,lower.tail=FALSE,log.p=TRUE))
  wi=pzi0+dnorm(zi0)*(hzi0-zi0)
  ne=sum(wi)
  W=diag(wi,nrow=n)
  Uno=matrix(1,nrow=n,ncol=1)
  Id=diag(1,n)
  ne.inv=Uno%*%t(Uno)/ne
  W.root=diag(sqrt(wi),nrow=n)
  Weig=W.root%*%(Id-W.root%*%ne.inv%*%W.root)%*%W.root
  Xt.Weig.X=t(X)%*%Weig%*%X
  
  if(det(Xt.Weig.X)<1e-10 || is.na(det(Xt.Weig.X))){
    X.sca=t(X)%*%(Id-Uno%*%t(Uno)/n)%*%X
    Xt.Weig.X.inv=n*solve(X.sca)
  }else{
    Xt.Weig.X.inv=ne*solve(Xt.Weig.X)
  }

  #For the inverse gamma:
  if(np>2) res=dmnorm(x=parms[-c(1,2)],mean=0,
                      varcov=exp(2*parms[1])*g*Xt.Weig.X.inv,log=TRUE)+ldinvgamma(g,0.5,0.5)
  #For g=1, comment the above line and uncomment this one:
  #if(np>2) res=dmnorm(x=parms[-c(1,2)],mean=0,
  #                    varcov=exp(2*parms[1])*1*Xt.Weig.X.inv,log=TRUE)
  return(res)
}


# Log-Normal with Conventional Prior, with g sim IGa(0.5,0.5)
log.post.lnormal.our=function (parms,g, Data,cens){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(log.prior.our(parms,g,cens,X=as.matrix(Data[,-(1:3)]))+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

# Log-Normal with Conventional Prior, with g sim IGa(0.5,0.5) in vectorial form, parms contains theta and g:
log.post.lnormal.our.vec=function (parms, Data,cens){
  p=length(parms)
  #g appears in the last position:
  z = (Data[,1] - cbind(Data[,-c(1:2)]) %*% parms[-c(1,p)])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(log.prior.our(parms[-p],parms[p],cens,X=as.matrix(Data[,-(1:3)]))+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

#Same posterior using the robust prior pi(g):
# Log-Normal with Conventional Prior, with g sim Pi(g) robust prior 
log.post.lnormal.robust=function (parms,g, Data,cens){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(log.prior.robust(parms,g,cens,X=as.matrix(Data[,-(1:3)]))+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

#Using as prior over g \sim pi(g) robust prior:
log.prior.robust=function(parms,g,cens,X){
  np=length(parms)
  n=length(cens)
  res=0
  zi0=(cens-parms[2])/exp(parms[1])
  pzi0=pnorm(zi0)
  Delta0=diag(pzi0,nrow=n)
  hzi0=exp(dnorm(zi0,log=TRUE)-pnorm(zi0,lower.tail=FALSE,log.p=TRUE))
  wi=pzi0+dnorm(zi0)*(hzi0-zi0)
  ne=sum(wi)
  W=diag(wi,nrow=n)
  Uno=matrix(1,nrow=n,ncol=1)
  Id=diag(1,n)
  ne.inv=Uno%*%t(Uno)/ne
  W.root=diag(sqrt(wi),nrow=n)
  Weig=W.root%*%(Id-W.root%*%ne.inv%*%W.root)%*%W.root
  Xt.Weig.X=t(X)%*%Weig%*%X
  
  if(det(Xt.Weig.X)<1e-10 || is.na(det(Xt.Weig.X))){
    X.sca=t(X)%*%(Id-Uno%*%t(Uno)/n)%*%X
    Xt.Weig.X.inv=n*solve(X.sca)
  }else{
    Xt.Weig.X.inv=ne*solve(Xt.Weig.X)
  }
  
  if(np>2) res=dmnorm(x=parms[-c(1,2)],mean=0,
                      varcov=exp(2*parms[1])*g*Xt.Weig.X.inv,log=TRUE)+lpi.g(g,n,k.gamma=np-2)
  return(res)
}



log.post.lnormal.our.gfixed=function (parms,g=1, Data,cens){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(log.prior.our.gfixed(parms,g=1,cens,X=as.matrix(Data[,-(1:3)]))+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}


log.post.lnormal.0.our=function (parms, Data){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}


##Function to calculate m1 given a set of covariantes
##Approximation marginal 1 under an alternative model:
#Prior for g: Robust prior.
marg1.allpara.robust=function(y,rel,ct,X,k.star,nsim.marg=1000,nsim.post=10000){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  n=length(y)
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=parms.hat1,fn=log.post.lnormal,control=list(fnscale=-1),Data=data1.aux,method=optimmethod,hessian=TRUE)
  parms.hat1=tt1$par
  var.prop=-solve(tt1$hessian)
  
  proposal=list(var=var.prop,scale=0.5)
  res=simul.post.robust.rw(log.post.lnormal.robust,proposal,start=parms.hat1,g.start=1,nsim.post,n=n,Data=data1.aux,cens=ct)
  res$par=res$par[-(1:(nsim.post*0.3)),]
  
  n.par=dim(res$par)[2]
  mu.prop=apply(res$par[,-n.par],2,mean)
  Sigma.prop=cov(res$par[,-n.par])
  
  thetas=rmvt(nsim.marg, delta = mu.prop, sigma =Sigma.prop,df=4)
  gs=r.robust.g(nsim.marg,n,k.gamma=ncov)
  
  cens=ct
  
  lnum=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x) log.post.lnormal.robust(thetas[x,],g=gs[x],Data=as.matrix(data1.aux),cens=cens))
  lnum=lnum+k.star 
  lden=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x){ dmvt(thetas[x,],delta=mu.prop, sigma =Sigma.prop,df=3,log=TRUE) + lpi.g(gs[x],n,k.gamma=ncov)})

  marg1.comp=mean(exp(lnum-lden),na.rm=TRUE)
  return(list(marg=marg1.comp,par.post=res$par))
}



#Calculation of the effective sample size:
calc.ne=function(parms,censtimes){
  zi0=(censtimes-parms[2])/exp(parms[1])
  pzi0=pnorm(zi0)
  hzi0=exp(dnorm(zi0,log=TRUE)-pnorm(zi0,lower.tail=FALSE,log.p=TRUE))
  wi=pzi0+dnorm(zi0)*(hzi0-zi0)
  return(sum(wi))
}

#Another alternative to calculate ne, not used in the paper:
calc.ne.alternative=function(parms,censtimes,X){
  zi0=(censtimes-X%*%parms[-1])/exp(parms[1])
  pzi0=pnorm(zi0)
  hzi0=exp(dnorm(zi0,log=TRUE)-pnorm(zi0,lower.tail=FALSE,log.p=TRUE))
  wi=pzi0+dnorm(zi0)*(hzi0-zi0)
  return(sum(wi))
}


#Function to simulate from the posterior distribution, it simulates from g, in this funciton using the inverse gamma employing MH and a RW:
simul.post.rw=function(logpost,proposal,start,g.start,m,...){
  pb = length(start)
  Mpar = array(0, c(m, pb+1))
  b = matrix(t(start))
  g=g.start
  lb = logpost(start,g.start, ...)
  lg=ldinvgamma(g.start,0.5,0.5)
  a = chol(proposal$var)
  scale = proposal$scale
  accept = 0
  for (i in 1:m) {
    bc = b + scale * t(a) %*% array(rnorm(pb), c(pb, 1))
    gc=rinvgamma(1,0.5,0.5)
    lbc = logpost(t(bc),gc, ...)
    lgc=ldinvgamma(gc,0.5,0.5)
    #The jumping probability:
    prob = exp(lbc - lb+lg-lgc)
    if (is.na(prob) == FALSE) {
      if (runif(1) < prob) {
        lb = lbc
        b = bc
        g=gc
        lg=lgc
        accept = accept + 1
      }
    }
    Mpar[i, ] = c(b,g)
  }
  accept = accept/m
  stuff = list(par = Mpar, accept = accept)
  return(stuff)
  
}

#Function to simulate from the posterior distribution when g is fixed, g=1
simul.post.rw.gfixed=function(logpost,proposal,start,g,m,...){
  pb = length(start)
  Mpar = array(0, c(m, pb))
  b = matrix(t(start))
  lb = logpost(start,g, ...)
  a = chol(proposal$var)
  scale = proposal$scale
  accept = 0
  for (i in 1:m) {
    bc = b + scale * t(a) %*% array(rnorm(pb), c(pb, 1))
    lbc = logpost(t(bc),g, ...)
    prob = exp(lbc - lb)
    if (is.na(prob) == FALSE) {
      if (runif(1) < prob) {
        lb = lbc
        b = bc
        accept = accept + 1
      }
    }
    Mpar[i, ] = c(b)
  }
  accept = accept/m
  stuff = list(par = Mpar, accept = accept)
  return(stuff)
  
}



#Function to simulate from the posterior distribution, it simulates from g, in this funciton using the robust prior employing MH and a RW:
simul.post.robust.rw=function(logpost,proposal,start,g.start,m,n,...){
  pb = length(start)
  Mpar = array(0, c(m, pb+1))
  b = matrix(t(start))
  g=g.start
  k.gamma=pb-2
  lb = logpost(start,g.start, ...)
  lg=lpi.g(g.start,n,k.gamma)
  a = chol(proposal$var)
  scale = proposal$scale
  accept = 0
  for (i in 1:m) {
    bc = b + scale * t(a) %*% array(rnorm(pb), c(pb, 1))
    gc=r.robust.g(1,n,k.gamma)
    lbc = logpost(t(bc),gc, ...)
    lgc=lpi.g(gc,n,k.gamma)
    #Jumping probability:
    prob = exp(lbc - lb+lg-lgc)
    if (is.na(prob) == FALSE) {
      if (runif(1) < prob) {
        lb = lbc
        b = bc
        g=gc
        lg=lgc
        accept = accept + 1
      }
    }
    Mpar[i, ] = c(b,g)
  }
  accept = accept/m
  stuff = list(par = Mpar, accept = accept)
  return(stuff)
  
}

#Simulation of g using Inverse Transform of Distribution:
r.robust.g=function(nsimul,n,k.gamma){
  unif=runif(nsimul,0,1)
  g.simul=(1+n)/(n*(k.gamma+1))*(1-unif)^(-2)-1/n
  return(g.simul)
}

#Density robust prior:
pi.g=function(g,n,k.gamma){
  if(g> (1+n)/(n*(k.gamma+1))-1/n) val.pi=0.5*sqrt((1+n)/((k.gamma+1)*n))*(g+1/n)^(-3/2)
  else val.pi=0
  return(val.pi)
}

#Logdensity robust prior:
lpi.g=function(g,n,k.gamma){
  if(g> (1+n)/(n*(k.gamma+1))-1/n) val.lpi=log(0.5)+0.5*(log(1+n)-log(k.gamma+1)-log(n))-3/2*log(g+1/n)
  else val.lpi=-Inf
  return(val.lpi)
}


#Posterior probability calculation, using a uniform prior over the model space (if scott=FALSE) or using the Scott-Berger prior if scott=TRUE
#using the BFs:             
prob.post=function(BF,mod.list,scott=TRUE){
  nmodels=length(mod.list)
  prob.post=rep(NA,nmodels+1)
  BF=c(1,BF)
  if(scott==FALSE){
    for(i in 1:(nmodels+1)){
      prob.post[i]=(1+sum(BF[-i]/BF[i]))^(-1)
    }
  }else{
    ncov=c(0,laply(mod.list,length))
    uncov=unique(ncov)
    p=max(ncov)
    redprior=1/(p+1)
    redprior=redprior/c(table(ncov))
    modprior=redprior[ncov+1]
    for(i in 1:(nmodels+1)){
      prob.post[i]=(1+sum(modprior[-i]*BF[-i]/(modprior[i]*BF[i])))^(-1)
    }
  }
  return(prob.post)
}

################################################
##Functions to calculate m1, when using as prior for beta a normal with var-covar defined as Sigma^A in the paper.
marg1.allpara.gprior.XtX=function(y,rel,cens,X,k.star,nsim.marg=1000,nsim.post=5000){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=parms.hat1,fn=log.post.lnormal,control=list(fnscale=-1),Data=data1.aux,method=optimmethod,hessian=TRUE)
  var.prop=-solve(tt1$hessian)
  
  parms.hat1=tt1$par
  
  proposal=list(var=var.prop,scale=1)
  #For g fixed, g=1
  res.1=simul.post.rw.gfixed(log.post.lnormal.XtX,proposal,start=parms.hat1,g=1,m=nsimul,Data=as.matrix(data1.aux),cens=cens)
  res.1$par=res.1$par[-(1:(nsimul*0.3)),]
  
  mu.prop=apply(res.1$par,2,mean)
  Sigma.prop=cov(res.1$par)
  
  thetas=rmvt(nsim.marg, delta = mu.prop, sigma =Sigma.prop,df=3)
  lnum=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x) log.post.lnormal.XtX(thetas[x,],g=1,Data=as.matrix(data1.aux),cens=cens))
  lnum=lnum+k.star 
  lden=apply(as.matrix(1:nsim.marg,ncol=nsim.marg),1,function(x){ dmvt(thetas[x,],delta=mu.prop, sigma =Sigma.prop,df=3,log=TRUE)})
  
  marg1.comp=mean(exp(lnum-lden),na.rm=TRUE)
  return(list(marg=marg1.comp,par.post=res.1$par))
}


# Log-Normal with Conventional Prior, and var-covar matrix Sigma^A
log.post.lnormal.XtX=function (parms,g, Data,cens){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(log.prior.XtX(parms,g,cens,X=as.matrix(Data[,-(1:3)]))+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}


# Log-prior with var-matrix Sigma^A:
log.prior.XtX=function(parms,g=1,cens,X){
  np=length(parms)
  n=length(cens)
  XtX=t(X)%*%X
  XtX.inv=n*solve(XtX)
  g=1
  
  if(np>2) res=dmnorm(x=parms[-c(1,2)],mean=0,
                      varcov=exp(2*parms[1])*g*XtX.inv,log=TRUE)
  return(res)
}


########################################
#Aprox using Laplace:
##Function to calculate m1 given a set of covariates
#Prior for g: Robust prior.

#Log-Normal with Conventional Prior, with g sim Pi(g) robust prior; with Sigma^M fixed for beta0 and sigma calculate under the null:
log.post.lnormal.robust.sigmaMfixed=function (g,parms,Data,cens,SigmaM){
  p=length(parms)
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-c(1)])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  V=g*SigmaM
  return(dmnorm(x=parms[-c(1,2)],mean=0,
                varcov=V,log=TRUE)+lpi.g(g,n,k.gamma=p-2)+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

post.lnormal.robust.sigmaMfixed.tointe.aux=function(g,parms,Data,cens,SigmaM,k.star=NULL){exp(log.post.lnormal.robust.sigmaMfixed(g=g,parms=parms,Data=Data,cens=cens,SigmaM=SigmaM)+k.star)}
post.lnormal.robust.sigmaMfixed.tointe=function(g,parms,Data,cens,SigmaM,k.star=NULL){ mapply(post.lnormal.robust.sigmaMfixed.tointe.aux,g,MoreArgs=list(parms=parms,Data=Data,cens=cens,SigmaM=SigmaM,k.star=k.star))}

#Log Posterior distribution for Log-Normal with Conventional Prior, where we integrate g with respecto to the g-prior
log.post.lnormal.robust.sigmaMfixed.g.inte=function(parms,Data,cens,SigmaM,k.star=NULL){
  n=length(cens)
  k.gamma=length(parms)-2
  tmp=integrate(post.lnormal.robust.sigmaMfixed.tointe,lower=(1+n)/(n*(k.gamma+1))-1/n,upper=Inf,parms=parms,Data=Data,cens=cens,SigmaM=SigmaM,k.star)
  return(log(tmp$value))
}

#Log-Normal with Inverse-Gamma Prior, with g sim InverseGamma(g); with Sigma^M fixed for beta0 and sigma calculate under the null:
log.post.lnormal.igamma.sigmaMfixed=function (parms,g,Data,cens,SigmaM){
  p=length(parms)
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-c(1)])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(dmnorm(x=parms[-c(1,2)],mean=0,
                varcov=g*SigmaM,log=TRUE)+ldinvgamma(g,0.5,0.5)+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

#Same function than before but written in vectorial form, all the parameters are in parms (also g).
log.post.lnormal.igamma.sigmaMfixed.vec=function (parms,Data,cens,SigmaM){
  p=length(parms)
  z = (Data[,1] - cbind(Data[,-c(1:2)]) %*% parms[-c(1,p)])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(dmnorm(x=parms[-c(1,2,p)],mean=0,
                varcov=parms[p]*SigmaM,log=TRUE)+ldinvgamma(parms[p],0.5,0.5)+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

############################################################
#With g fixed to 1
log.post.lnormal.gfixed.sigmaMfixed=function (parms,g=1, Data,cens,SigmaM){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(dmnorm(x=parms[-c(1,2)],mean=0,
                varcov=SigmaM,log=TRUE)+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

############################################################
#With g fixed to 1, written in vectorial form, all the parameters are in parms (also g).
log.post.lnormal.gfixed.sigmaMfixed.vec=function (parms,Data,cens,SigmaM){
  z = (Data[,1] - cbind(Data[,-c(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(dmnorm(x=parms[-c(1,2)],mean=0,
                varcov=SigmaM,log=TRUE)+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

#Function to calculate Sigma^M
calculo.SigmaM=function(parms,cens,X){
  n=length(cens)
  res=0
  zi0=(cens-parms[2])/exp(parms[1])
  pzi0=pnorm(zi0)
  Delta0=diag(pzi0,nrow=n)
  hzi0=exp(dnorm(zi0,log=TRUE)-pnorm(zi0,lower.tail=FALSE,log.p=TRUE))
  wi=pzi0+dnorm(zi0)*(hzi0-zi0)
  ne=sum(wi)
  W=diag(wi,nrow=n)
  Uno=matrix(1,nrow=n,ncol=1)
  Id=diag(1,n)
  ne.inv=Uno%*%t(Uno)/ne
  W.root=diag(sqrt(wi),nrow=n)
  Weig=W.root%*%(Id-W.root%*%ne.inv%*%W.root)%*%W.root
  Xt.Weig.X=t(X)%*%Weig%*%X
  
  if(det(Xt.Weig.X)<1e-10 || is.na(det(Xt.Weig.X))){
    X.sca=t(X)%*%(Id-Uno%*%t(Uno)/n)%*%X
    Xt.Weig.X.inv=n*solve(X.sca)
  }else{
    Xt.Weig.X.inv=ne*solve(Xt.Weig.X)
  }
  SigmaM=exp(2*parms[1])*Xt.Weig.X.inv
  return(SigmaM)
}

##Marginal m1 calculation, using the robust prior and the Laplace approximation with fixed SigmaM  (for beta0 and sigma calculate under the null):
marg1.allpara.robust.laplace=function(y,rel,ct,X,k.star=NULL,nsim.marg=NULL,nsim.post=NULL,beta0.hat,lsigma.hat,g0.init=0){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  n=length(y)
  SigmaM.hat=calculo.SigmaM(parms=c(lsigma.hat,beta0.hat),cens=ct,X)
  
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=parms.hat1,fn=log.post.lnormal.robust.sigmaMfixed.g.inte,Data=data1.aux,
            SigmaM=SigmaM.hat,cens=ct,k.star=k.star,hessian=TRUE,method="Nelder-Mead",control=list(fnscale=-1))
  hessi=-tt1$hessian
  
  sd.par=sqrt(diag(solve(hessi)))
   
  map.coef=as.vector(tt1$par)
  p=dim(hessi)[1]
  lm1.lapl=log.post.lnormal.robust.sigmaMfixed.g.inte(parms=map.coef,Data=data1.aux,cens=ct,SigmaM=SigmaM.hat,k.star=k.star)+p/2*log(2*pi)-0.5*log(det(hessi))
  return(list(marg=lm1.lapl,par.post=NULL))
}


##Marginal m1 calculation, using the inverse gamma prior and the Laplace approximation with fixed SigmaM  (for beta0 and sigma calculate under the null):
marg1.allpara.igamma.laplace=function(y,rel,ct,X,k.star=NULL,nsim.marg=NULL,nsim.post=NULL,beta0.hat,lsigma.hat){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  n=length(y)
  SigmaM.hat=calculo.SigmaM(parms=c(lsigma.hat,beta0.hat),cens=ct,X)
  
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=c(parms.hat1,1),fn=log.post.lnormal.igamma.sigmaMfixed.vec,control=list(fnscale=-1),Data=data1.aux,SigmaM=SigmaM.hat,cens=ct,hessian=TRUE)
  
  hessi=-tt1$hessian
  
  map.coef=as.vector(tt1$par)
  p=dim(hessi)[1]
 
  lm1.lapl=log.post.lnormal.igamma.sigmaMfixed(parms=map.coef[-p],g=map.coef[p],Data=data1.aux,cens=ct,SigmaM=SigmaM.hat)+p/2*log(2*pi)+log(det(hessi)^(-0.5))
  return(list(marg=lm1.lapl,par.post=NULL))
}


##Marginal m1 calculation, using the inverse gamma prior and the Laplace approximation without fixing SigmaM:
marg1.allpara.igamma.laplace2=function(y,rel,ct,X,k.star=NULL,nsim.marg=NULL,nsim.post=NULL,beta0.hat,lsigma.hat){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  n=length(y)
  
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=c(parms.hat1,1),fn=log.post.lnormal.our.vec,control=list(fnscale=-1),Data=data1.aux,cens=ct,hessian=TRUE)
  
  hessi=-tt1$hessian
  
  map.coef=as.vector(tt1$par)
  p=dim(hessi)[1]
  
  lm1.lapl=log.post.lnormal.our(parms=map.coef[-p],g=map.coef[p],Data=data1.aux,cens=ct)+p/2*log(2*pi)+log(det(hessi)^(-0.5))
  return(list(marg=lm1.lapl,par.post=NULL))
}

##Marginal m1 calculation, with g=1 fixed and the Laplace approximation  with fixed SigmaM  (for beta0 and sigma calculate under the null):
marg1.allpara.gfixed.laplace=function(y,rel,ct,X,k.star=NULL,nsim.marg=NULL,nsim.post=NULL,beta0.hat,lsigma.hat){
  ncov=ncol(X)
  data1.aux=as.matrix(cbind(y,rel,1,X))
  n=length(y)
  SigmaM.hat=calculo.SigmaM(parms=c(lsigma.hat,beta0.hat),cens=ct,X)
  
  colnames(data1.aux)=c("y","rel","Intercp",colnames(X))
  tt=survreg(Surv(exp(y),rel)~data1.aux[,-c(1:3)],dist="lognormal")
  parms.hat1=c(log(tt$scale),tt$coeff)
  
  tt1=optim(par=parms.hat1,fn=log.post.lnormal.gfixed.sigmaMfixed.vec,control=list(fnscale=-1),Data=data1.aux,SigmaM=SigmaM.hat,cens=ct,hessian=TRUE)
  
  hessi=-tt1$hessian
  
  map.coef=as.vector(tt1$par)
  p=dim(hessi)[1]
  
  lm1.lapl=log.post.lnormal.gfixed.sigmaMfixed(parms=map.coef,Data=data1.aux,cens=ct,SigmaM=SigmaM.hat)+p/2*log(2*pi)+log(det(hessi)^(-0.5))
  return(list(marg=lm1.lapl,par.post=NULL))
}
