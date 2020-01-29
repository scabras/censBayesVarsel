## install.packages(c("igraph0","BayesTree", "ggplot2","Hmisc","MCMCpack","LearnBayes","evd","plyr","dummies","bayespack"))
# install.packages(c("sp","numDeriv", "fields","rgl","mvtnorm","multicore","pixmap"))

library(MCMCpack)
library(LearnBayes)
library(evd)
library(plyr)
library(compiler)
library(parallel)
library(dummies)
library(survival)
#library(bayespack)

### Constants
l2p=log(2)+log(pi)
optimmethod=c("Nelder-Mead","CG")[1]

## Log Likelihoods and Kernel Posteriors with:
# b fraction (b=1 for usual)
# prior: unused

##Simulations from the generalized gumbel with k=2
rggumbel=function(nrep,nosirve1=0,nosirve2=1){
  ts=rgamma(nrep,shape=2,rate=1)
  return(ws=log(ts))
}
qggumbel=function(pcens,nosirve1=0,nosirve2=1) return(log(qgamma(pcens,shape=2,rate=1)))


# Log-Generalized Gamma k=2
###################################
log.post.ggamma.nc=function (parms, Data,b,prior=1){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]+2*z-exp(z)
  log.S=-exp(z)+log(1+exp(z))
  return(sum(Data[,2] * log.f*b[1] + (1 - Data[,2]) * log.S*b[2]))
}

log.prior.cp.gprior=function(parms,nunc,XtX.inv){
  np=length(parms)
  res=0
  if(np>2) res=dmnorm(x=parms[-c(1,2)],mean=0,
                      varcov=exp(2*parms[1])*nunc*XtX.inv,log=TRUE)
  return(res)
}
log.prior.cp.zs=function(parms,nunc,XtX.inv){
  np=length(parms)
  res=0
  if(np>2) res=dmt(x=parms[-c(1,2)],
                      S=exp(2*parms[1])*nunc*XtX.inv,df=1,log=TRUE)
  return(res)
}


# Log-Normal with Conventional Prior
log.post.lnormal.cp.nc=function (parms, Data,b,prior=1,XtX.inv=NULL,nunc=NULL){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(log.prior.cp(parms,nunc,XtX.inv)+sum(Data[,2] * log.f + (1 - Data[,2]) * log.S))
}

# Log-Normal
log.post.lnormal.nc=function (parms, Data,b,prior=1,XtX.inv=NULL,nunc=NULL){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1]-0.5*z^2
  log.S=pnorm(z,lower.tail=FALSE,log.p=TRUE)
  return(sum(Data[,2] * log.f*b[1] + (1 - Data[,2]) * log.S*b[2]))
}
# Weibull
log.post.weibull.nc=function (parms, Data,b,prior=1,XtX.inv=NULL,nunc=NULL){
  z = (Data[,1] - cbind(Data[,-(1:2)]) %*% parms[-1])/exp(parms[1])
  log.f = -parms[1] +(z - exp(z))
  log.S = -exp(z)
  return(sum(Data[,2] * log.f*b[1] + (1 - Data[,2]) * log.S*b[2]))
}

log.post.weibull=cmpfun(log.post.weibull.nc)
log.post.lnormal=cmpfun(log.post.lnormal.nc)
log.post.lnormal.cp=cmpfun(log.post.lnormal.cp.nc)
log.post.ggamma=cmpfun(log.post.ggamma.nc)


############# Posterior Calculation using MCMC
post=function(log.post,nsimul,parms.init,b,dat,nunc,XtX.inv){
  fit=optim(parms.init,log.post,b=b,Data=dat,XtX.inv=XtX.inv,nunc=nunc,control=list(fnscale=-1),method=optimmethod,hessian=TRUE)
  var.prop=-solve(fit$hessian)
  proposal=list(var=var.prop,scale=2)
  bayesfit=rwmetrop(log.post,proposal,fit$par,nsimul,b=b,Data=dat,XtX.inv=XtX.inv,nunc=nunc)
  stuff=list(var.prop=(proposal$scale^2)*proposal$var,par=bayesfit$par)
  return(stuff)
}



# Calculation of log Marginal di Chibb-Jeliazkov
log.marg.nc=function (log.post, parms.hat, par, var.prop, dat, b,nunc,XtX.inv){
  nsimul=nrow(par)
  #calcolo della posterior approssimata
  lbc = log.post(parms.hat, dat, b,XtX.inv=XtX.inv,nunc=nunc)
  lb=aaply(par,1,log.post,Data=dat,b=b,XtX.inv=XtX.inv,nunc=nunc)
  num=mean(exp(lbc-lb+dmnorm(par,parms.hat,var.prop,log=T)))
  parms.j=mvrnorm(nsimul,mu=parms.hat,Sigma=var.prop)
  lbc.den = aaply(parms.j,1,log.post,Data=dat,b=b,XtX.inv=XtX.inv,nunc=nunc)
  lb.den = lbc	 
  den=mean(exp(lbc.den-lb.den))
  return(log.post(parms.hat,Data=dat,b=b,XtX.inv=XtX.inv,nunc=nunc)-log(num)+log(den))
}
log.marg=cmpfun(log.marg.nc)

# Calculation of log Marginal Laplace using banint of package BayesPack
log.marg.laplace.bp=function (log.post, parms.hat, dat, b,XtX.inv,nunc){
  mlp=function(parms) log.post(parms,Data=dat,b=b,XtX.inv=XtX.inv,nunc=nunc)
  res=banint(mlp,mode=parms.hat,optimMethod = optimmethod,method="Monte Carlo")$nrmcon[1]
  return(res)
}

# Calculation of log Marginal Laplace
log.marg.laplace.nc=function (log.post, parms.hat, dat, b,XtX.inv,nunc){
  d=length(parms.hat)
  N=nrow(dat)
  fit=optim(parms.hat,log.post,b=b,Data=dat,XtX.inv=XtX.inv,nunc=nunc,control=list(fnscale=-1),method=optimmethod,hessian=TRUE)
  M=fit$value+l2p*d*0.5-log(abs(det(fit$hessian)))*0.5
  return(M)
}

log.marg.laplace=log.marg.laplace.nc

nunc.func.nu=function(x) sum(x)
nunc.func.n=function(x)  length(x)
nunc.func.uno=function(x)  1

func.XtX.inv.n=function(rel,XM,y=NULL,model=NULL) solve(t(XM)%*%XM)
func.XtX.inv.nu=function(rel,XM,y=NULL,model=NULL) solve(t(XM[rel==1,])%*%XM[rel==1,])

func.XtX.obsfisher.old=function(rel,XM,y=NULL){
  tt=survreg(Surv(exp(y),rel)~ XM,dist="lognormal")
  SIG=tt$var[-c(1,nrow(tt$var)),-c(1,nrow(tt$var))]/tt$scale^2
  return(SIG)
} 

func.XtX.obsfisher=function(rel,XM,y=NULL){
  tt=survreg(Surv(exp(y),rel)~ XM,dist="lognormal")
  ic=rel==0
  iu=!ic
  XMw.u=iu*XM
  XtX.uw=t(XMw.u)%*%XMw.u
  resid=c((y-tt$coef%*%t(cbind(rep(1,n),XM)))/tt$scale)
  ww=sqrt(2*pi)*(1-pnorm(resid))
  ww=1/ww*(exp(-resid^2)/ww-resid*exp(-resid^2/2))*ic
  XMw.c=sqrt(ww)*XM
  XtX.cw=t(XMw.c)%*%XMw.c
  I.F=XtX.uw/sum(iu)+XtX.cw/sum(ic)
  SIG=solve(I.F)
  return(SIG)
}



############# Log of BF 1 against 0 based on data0 and data1
log.bf.parziale=function(data0,data1,nsimul,parms.init.used,b,log.post,use.MCMC=NULL){
  parms.hat0=parms.init.used[[1]]
  parms.hat1=parms.init.used[[2]]
  
  nunc=nunc.func(data1[,2])
  XtX.inv=func.XtX.inv(rel=data1[,2],XM=cbind(data1[,-(1:3)]),y=data1[,1])
    
    
  
  if(use.MCMC){ ## Use MCMC
    res=post(log.post,nsimul,parms.hat0,dat=data0,b=b,XtX.inv=NULL,nunc=NULL)
    res$par=res$par[-(1:(nsimul*0.3)),]
    parms.hat0=aaply(res$par,2,median)
    log.marg0=log.marg(log.post,parms.hat0,res$par,res$var.prop,dat=data0,b=b)
    res=post(log.post,nsimul,parms.hat1,dat=data1,b=b,nunc=nunc,XtX.inv=XtX.inv)
    res$par=res$par[-(1:(nsimul*0.3)),]
    parms.hat1=aaply(res$par,2,median)
    log.marg1=log.marg(log.post,parms.hat1,res$par,res$var.prop,
                       dat=data1,b=b,XtX.inv=XtX.inv,nunc=nunc)
  }
  #log.marg1-log.marg0<0
  if(use.MCMC==FALSE){ ## Use Laplace
    log.marg0=log.marg.laplace(log.post=log.post,
                               parms.hat=parms.hat0,dat=data0,b=b,XtX.inv=NULL,nunc=NULL)
    log.marg1=log.marg.laplace(log.post=log.post,
                               parms.hat=parms.hat1,dat=data1,b=b,XtX.inv=XtX.inv,nunc=nunc)
  }
  return(log.marg1-log.marg0)
}

############# Data Simulation according to the model specified
sim.dati=function(n,beta.true,p.cens = 0.3,sigma,mu,model="lnorm"){
  parms.true = c(log(sigma), mu, beta.true)
  p = length(beta.true)
  ones = rep(1,n)
  x = cbind(ones, array(rnorm(n*p,sd=3), dim=c(n,p)))
  
  #if(p==1) x[,2]=x[,2]-mean(x[,2]) else x[,-1]=apply(x[,-1],2,function(xx) xx-mean(xx))
  mu1=x%*%parms.true[-1]
  if(model=="lnorm"){
    rmodel=rnorm
    mu2 = mu1-sqrt(2*sigma^2)*qnorm(p.cens,0,1)
  }
  if(model=="weibull"){
    rmodel=rgumbel
    
    lambda = exp(-(mu1)/exp(parms.true[1]))
    alpha = 1/exp(parms.true[1])
    rate.weib = lambda*p.cens/(1-p.cens)
  }  
  if(model=="ggamma"){
    rmodel=rggumbel
    mu2 = mu1-qggumbel(p.cens,0,1)
    }
  if((model=="weibull")|(model=="lnorm")) y = mu1-exp(parms.true[1])*rmodel(n,0,1)
  if(model=="ggamma")   y = mu1+exp(parms.true[1])*rmodel(n,0,1)
  ## Cens Part
  ncens = round(n*p.cens,0)
  relapse = ones
  ids = rep(FALSE,n)
  while(sum(ids)<ncens){
    if(model=="weibull") censtimes = log(rweibull(n, shape=alpha, scale=(1/rate.weib)^(1/alpha)))
    if(model=="lnorm") censtimes = rnorm(n,mu2,sigma)
    if(model=="ggamma") censtimes = mu2+rmodel(n,0,1)
    ids = censtimes<y
  }
  if(sum(ids)>0){
    #ind.cens = sample(which(ids), ncens)
    #y[ind.cens] = censtimes[ind.cens]
    y[ids]=censtimes[ids]
    relapse[ids] = 0
  }
  dat=cbind(y,relapse,x)
  parms.init = list(init0=parms.true[c(1,2)], init1=parms.true)
  return(list(dat=dat,parms.init=parms.init,censtimes=censtimes))
}

#sim.dati.detcens(n,beta.true,p.cens,sigma=sigma.true,mu=mu.true,model=model)
############# Data Simulation according to the model specified
sim.dati.detcens=function(n,beta.true,p.cens = 0.3,sigma,mu,model="lnorm"){
  parms.true = c(log(sigma), mu, beta.true)
  p = length(beta.true)
  ones = rep(1,n)
  x = cbind(ones, array(rnorm(n*p,sd=3), dim=c(n,p)))
  
  if(p==1) x[,2]=x[,2]-mean(x[,2]) else x[,-1]=apply(x[,-1],2,function(xx) xx-mean(xx))
  mu1=x%*%parms.true[-1]
  
  ind.ord=sort(mu1,index.return=TRUE)$ix
  mu1=mu1[ind.ord]
  x=x[ind.ord,]
  
  ncens = round(n*p.cens,0)
  id.cens=c(1:(ncens/2),floor((n-ncens/2+1)):n)

  p.cens2=0.99
  p.cens3=0.01
  mu2=rep(NA,length(mu1))
  if(model=="lnorm"){
     rmodel=rnorm
     mu2[id.cens]=mu1[id.cens]-sqrt(2*sigma^2)*qnorm(p.cens2,0,1)
     mu2[-id.cens] = mu1[-id.cens]-sqrt(2*sigma^2)*qnorm(p.cens3,0,1)
   }
   if(model=="weibull"){
     rmodel=rgumbel
     
     lambda = exp(-(mu1)/exp(parms.true[1]))
     alpha = 1/exp(parms.true[1])
     rate.weib = lambda*p.cens2/(1-p.cens2)
   }  
   if(model=="ggamma"){
     rmodel=rggumbel
     mu2 = mu1-qggumbel(p.cens2,0,1)
   }
   if((model=="weibull")|(model=="lnorm")) y = mu1-exp(parms.true[1])*rmodel(n,0,1)
   if(model=="ggamma")   y = mu1+exp(parms.true[1])*rmodel(n,0,1)
   ## Cens Part
   relapse = ones
   ids = rep(FALSE,n)
   if(model=="weibull") censtimes = log(rweibull(n, shape=alpha, scale=(1/rate.weib)^(1/alpha)))
   if(model=="lnorm") censtimes = rnorm(n,mu2,sigma)
   if(model=="ggamma") censtimes = mu2+rmodel(n,0,1)
   ids = censtimes<y
   relapse[ids] = 0
   y[ids] = censtimes[ids]
    dat=cbind(y,relapse,x)
   parms.init = list(init0=parms.true[c(1,2)], init1=parms.true)
   return(list(dat=dat,parms.init=parms.init,censtimes=censtimes))
}

############# Data Simulation according to the model specified
sim.dati.detcens.mahala=function(n,beta.true,p.cens = 0.3,sigma,mu,model="lnorm"){
  parms.true = c(log(sigma), mu, beta.true)
  p = length(beta.true)
  ones = rep(1,n)
  x = cbind(ones, array(rnorm(n*p,sd=3), dim=c(n,p)))
  
  if(p==1){ x[,2]=x[,2]-mean(x[,2])} else{x[,-1]=apply(x[,-1],2,function(xx) xx-mean(xx))}
  
  
  if(sum(beta.true)>0){
    covari=as.matrix(x[,-1])
    var.x=var(covari[,which(beta.true!=0)])
    maha.dist=mahalanobis(as.matrix(covari[,which(beta.true!=0)],nrow=n), rep(0,sum(beta.true)), cov=var.x)
    ind.ord=sort(maha.dist,index.return=TRUE)$ix}else{ind.ord=1:n}
  
  mu1=x%*%parms.true[-1]
  mu1=mu1[ind.ord]
  x=x[ind.ord,]
  
  ncens = round(n*p.cens,0)
  id.cens=(n-ncens+1):n
  
  if(model=="lnorm"){
    rmodel=rnorm
    mu2 = mu1-sqrt(2*sigma^2)*qnorm(p.cens,0,1)
    
   }
  if(model=="weibull"){
    rmodel=rgumbel
    
    lambda = exp(-(mu1)/exp(parms.true[1]))
    alpha = 1/exp(parms.true[1])
    rate.weib = lambda*p.cens2/(1-p.cens2)
  }  
  if(model=="ggamma"){
    rmodel=rggumbel
    mu2 = mu1-qggumbel(p.cens2,0,1)
  }
  if((model=="weibull")|(model=="lnorm")) y = mu1-exp(parms.true[1])*rmodel(n,0,1)
  if(model=="ggamma")   y = mu1+exp(parms.true[1])*rmodel(n,0,1)
  ## Cens Part
  relapse = ones
  ids = rep(FALSE,n)
  if(model=="weibull") censtimes = log(rweibull(n, shape=alpha, scale=(1/rate.weib)^(1/alpha)))
  if(model=="lnorm") censtimes = rnorm(n,mu2,sigma)
  if(model=="ggamma") censtimes = mu2+rmodel(n,0,1)
  #ids = censtimes<y
  ids=id.cens
  relapse[ids] = 0
  y[ids] = censtimes[ids]
  dat=cbind(y,relapse,x)
  parms.init = list(init0=parms.true[c(1,2)], init1=parms.true)
  return(list(dat=dat,parms.init=parms.init,censtimes=censtimes))
}

############# Data Simulation according to the model specified
sim.dati.known.cens=function(n,beta.true,p.cens = 0.3,sigma,mu,model="lnorm"){
  parms.true = c(log(sigma), mu, beta.true)
  p = length(beta.true)
  ones = rep(1,n)
  x = cbind(ones, array(rnorm(n*p,sd=3), dim=c(n,p)))
  
  if(p==1) x[,2]=x[,2]-mean(x[,2]) else x[,-1]=apply(x[,-1],2,function(xx) xx-mean(xx))
  mu1=x%*%parms.true[-1]
  if(model=="lnorm"){
    rmodel=rnorm
    mu2 = mu1-sqrt(2*sigma^2)*qnorm(p.cens,0,1)
  }
  if(model=="weibull"){
    rmodel=rgumbel
    
    lambda = exp(-(mu1)/exp(parms.true[1]))
    alpha = 1/exp(parms.true[1])
    rate.weib = lambda*p.cens/(1-p.cens)
  }  
  if(model=="ggamma"){
    rmodel=rggumbel
    mu2 = mu1-qggumbel(p.cens,0,1)
  }
  if((model=="weibull")|(model=="lnorm")) y = mu1-exp(parms.true[1])*rmodel(n,0,1)
  if(model=="ggamma")   y = mu1+exp(parms.true[1])*rmodel(n,0,1)
  ## Cens Part
  ncens = round(n*p.cens,0)
  relapse = ones
  ids = rep(FALSE,n)
  
   if(model=="weibull") censtimes = log(rweibull(n, shape=alpha, scale=(1/rate.weib)^(1/alpha)))
  if(model=="lnorm") censtimes = rnorm(n,mu2,sigma)
  if(model=="ggamma") censtimes = mu2+rmodel(n,0,1)
  ids = censtimes<y
  
  y[ids] = censtimes[ids]
  relapse[ids] = 0
  
  dat=cbind(y,relapse,x)
  parms.init = list(init0=parms.true[c(1,2)], init1=parms.true)
  return(list(dat=dat,parms.init=parms.init,censtimes=censtimes))
}





############# Calculate the model indices based on only additive effects of the covariates
index.models=function(numvar,onlyone=FALSE){  
  indici = seq(1,numvar)
  ind.betas = NULL
  for(j in 1:(numvar*(1-onlyone)+1*(onlyone))){
    comb.beta = combn(indici, j)
    for(i in 1:dim(comb.beta)[2]) ind.betas=c(ind.betas,list(comb.beta[,i]))
  }
  return(ind.betas)
}

############# Probability distribution of Nt
prob.nstar=function(nstar,n,ncens,p) exp(lchoose(ncens,nstar-p)+lchoose(n-ncens-1,p-1)+lfactorial(nstar-1)+log(n-ncens)-(lfactorial(n)-lfactorial(n-nstar)))

# Probability distribution of B:
# use="all" : all B and probabilities;
# use="mode" : mode of B
# use="median" : median of B
prob.b=function(n,ncens,p,use,toplot=FALSE){
  nstar.possible=seq(p,ncens+p)
  nstar.dens=aaply(nstar.possible,1,prob.nstar,n=n,ncens=ncens,p=p)
  if(toplot) barplot(height=nstar.dens,names.arg=nstar.possible,xlab="nt",ylab="de",main="")
  if(use=="all"){
    b.frac=cbind(p/(n-ncens),(nstar.possible-p)/ncens)
    b.frac.prob=nstar.dens}
  if(use=="mode"){
    nstar.mode=nstar.possible[nstar.dens==max(nstar.dens)][1]
    b.frac.prob=1
    b.frac=cbind(p/(n-ncens),(nstar.mode-p)/ncens)}
if(use=="median"){
    nstar.median=min(nstar.possible[cumsum(nstar.dens)>=0.5])[1]
    b.frac.prob=1
    b.frac=cbind(p/(n-ncens),(nstar.median-p)/ncens)}
  return(list(B=b.frac,B.prob=b.frac.prob))
}


############# Return the indices for a Training sample
training.sample.conv.prior=function(dati,num.par=2){
  n=nrow(dati)
  cens=dati[,2]
  #spazio campionario
  camp=1:n
  ind.train=NULL
  while(sum(cens[ind.train])<num.par){
    ind.train=c(ind.train,sample(camp,size=1))
    camp=camp[!(camp%in%ind.train)]
  }
  return(ind.train)
}




######## Bayes Factors ########
# Fractional Bayes Factor with specified fractions Bs and probabiities
# for a certain model specified in ind.beta
fbf.ind.beta=function(data0,data1,Bs,Bmode,use.pericchi2005.fbf=FALSE,
                          parms.init.used,nsimul,log.post,force.L=NULL,use.MCMC=NULL,log.bf.all=NULL){
  B=Bs$B
  B.prob=Bs$B.prob
  
  if(use.pericchi2005.fbf){
    tmp=aaply(B,1,function(x){
      log.marg0=log.post(parms=parms.init.used[[1]],Data=data0,b=1-x)
      log.marg1=log.post(parms=parms.init.used[[2]],Data=data1,b=1-x)
      m1=length(parms.init.used[[2]])
      m0=length(parms.init.used[[1]])
      return( (log.marg1-log.marg0)+0.5*(m1-m0)*log(x) )})
  }else{
    tmp=aaply(B,1,function(x) log.bf.parziale(data0,data1,nsimul,parms.init.used=parms.init.used,
                                              b=x,log.post=log.post,use.MCMC=use.MCMC))
  }
  a=(1-use.pericchi2005.fbf)
  b=c(-1,1)[1+use.pericchi2005.fbf]
  resp=exp(log.bf.all*a+b*tmp)
  nvalid=(!(abs(resp)<Inf))|is.na(resp)
  resp[nvalid]=sign(resp[nvalid])*2*max(resp[!nvalid],na.rm=TRUE)
  mfbf=log(sum(resp*B.prob))
  whichmode=B.prob==max(B.prob)
  fbfmode=log(resp[whichmode])
  names(mfbf)=names(fbfmode)=NULL
  return(list(mfbf=mfbf,fbfmode=fbfmode))
}

fbf.ind.betalnorm=function(data0,data1,Bs,Bmode,use.pericchi2005.fbf=FALSE,
                           parms.init.used,nsimul,log.post,force.L=NULL,use.MCMC=NULL,log.bf.all=NULL){
  B=Bs$B
  B.prob=Bs$B.prob
  tmp=aaply(B,1,function(x){
    Ri=lm(y~1,data=data.frame(data0))
    Ri=sum(Ri$residuals^2)
    Rj=lm(y~.,data=data.frame(data1[,-2]))
    Rj=sum(Rj$residuals^2)
    b=x
    kj=ncol(data1)-2
    ki=ncol(data0)-2
    qj=0
    qi=0
    r=b*nrow(data0)
    lbfb=(lgamma(0.5*(r-ki))+
            lgamma(0.5*(n-kj))-lgamma(0.5*(r-kj))-
            lgamma(0.5*(n-ki))+(n-r)*0.5*(log(Ri)-log(Rj)))
    
    return(lbfb)})
  mfbf=log(sum(exp(tmp)*B.prob))
  fbfmode=tmp[B==Bmode$B]
  names(mfbf)=names(fbfmode)=NULL
  return(list(mfbf=mfbf,fbfmode=fbfmode))
}

# Intrindic Bayes Factor with specified fractions Bs and probabiities
# for a certain model specified in ind.beta
# force.L : if not null force to use exactly L training sample 
ibf.ind.beta=function(data0,data1,Bs,parms.init.used,nsimul,log.post,use.pericchi2005.fbf=NULL,
                      force.L=NULL,use.MCMC=NULL,Bmode=NULL,log.bf.all=NULL){
  if(is.null(force.L)){L=round(Bs[[1]]$B[1]*nrow(data0)^2,0)}else{L=force.L}
  tmp=aaply(1:L,1,function(x){
    lbfl="err"
    while(is.character(lbfl)){
      ts=training.sample(data1)
      lbfl=try(log.bf.parziale(data0[ts,],data1[ts,],nsimul,parms.init.used,
                               b=c(1,1),log.post,use.MCMC=use.MCMC),silent=TRUE)
    }
    return(lbfl)})
  mibf=log.bf.all-median(tmp,na.rm=TRUE)
  etmp=exp(-tmp)
  etmp=etmp[etmp<Inf]
  aibf=log.bf.all+log(mean(etmp,na.rm=TRUE))
  names(mibf)=names(aibf)=NULL
  return(list(mibf=mibf,aibf=aibf))
}

# Return the indices for a Kaplan-Meieir Training sample
mts.kaplan.meier = function(mydata,p)
{
  # Data set under the null
  mydata = mydata[order(mydata[,1]),] # Sort data (same order as surv$time)
  # KM Estimator
  mysurv = survfit(Surv(exp(mydata[,1]), mydata[,2]) ~ 1,type="kaplan-meier") 
  n = nrow(mydata) # Sample size
  mts.size = p # MTS Sample Size
  probs = 1-mysurv$surv # CDF
  probs[n] = 1 # Adjusted CDF at the maximum
  probs = c(probs[1],probs[-1]-probs[-n]) # Point mass at each Observations
  not.censored.obs = which(mydata[,2]==1) # Not censored Obs
  elle = sample(x=not.censored.obs,size=mts.size,prob=probs[not.censored.obs]) 
  # The elle
  return(elle)
}

lgamma05=lgamma(0.5)

log.bf.parzialekm=function(data0,data1){
  k1 = dim(data1)[2]-2
  k0 = dim(data0)[2]-2
  r = k1-k0
  n = dim(data0)[1]
  t.med = sum(data0[,1])/n
  X1 = data1[,-c(1,2)]
  X0 = data0[,-c(1,2)]
  inv = solve(t(X1)%*%X1)
  beta1.hat = (inv%*%t(X1)%*%data1[,1])
  R0=data0[,1]-t.med
  R0 = c( t(R0)%*%(R0))
  R1=data1[,1]-X1%*%beta1.hat
  R1 = c( t(R1) %*% (R1))
  # BF i contro 0
  ibf.N.l = k1*lgamma05+lgamma((n-k1)*0.5) +0.5*log(det(t(X0)%*%X0))+((n-k0)*0.5)*log(R0)-(k0*lgamma05+lgamma((n-k0)*0.5)+0.5*log(det(t(X1)%*%X1))+((n-k1)*0.5)*log(R1) )
  return(ibf.N.l)
}

# Intrindic Bayes Factor KM (normal only) with specified fractions Bs and probabiities
# for a certain model specified in ind.beta
# force.L : if not null force to use exactly L training sample 
ibfkm.ind.beta=function(data0,data1,Bs,parms.init.used,nsimul,log.post,use.pericchi2005.fbf=NULL,
                      force.L=NULL,use.MCMC=NULL,Bmode=NULL,log.bf.all=NULL){
  if(is.null(force.L)){L=round(Bs[[1]]$B[1]*nrow(data0)^2,0)}else{L=force.L}
  tmp=aaply(1:L,1,function(x){
    lbfl="err"
    while(is.character(lbfl)){
      ts=mts.kaplan.meier(data0,p=ncol(data1)-1)
      lbfl=try(log.bf.parzialekm(data0[ts,],data1[ts,]),silent=TRUE)
    }
    return(lbfl)})
  mibf=log.bf.all-median(tmp,na.rm=TRUE)
  etmp=exp(-tmp)
  etmp=etmp[etmp<Inf]
  aibf=log.bf.all+log(mean(etmp,na.rm=TRUE))
  names(mibf)=names(aibf)=NULL
  return(list(mibf=mibf,aibf=aibf))
}

# BIC with specified fractions Bs and probabiities
# for a certain model specified in ind.beta
bic.ind.beta=function(data0,data1,Bs,parms.init.used,nsimul,log.post,force.L=NULL,use.pericchi2005.fbf=NULL,
                      use.MCMC=NULL,Bmode=NULL,log.bf.all=NULL){
  n.unc=sum(data0[,2])
  log0=log.post(parms=parms.init.used[[1]],Data=data0,b=c(1,1))
  log1=log.post(parms=parms.init.used[[2]],Data=data1,b=c(1,1))
  m0=length(parms.init.used[[1]])
  m1=length(parms.init.used[[2]])
  resbic=log1-log0+0.5*(m0-m1)*log(n.unc)
  names(resbic)=NULL
  return(resbic)
}


mppm.func=function(probs,ind.betas,nmodels=NULL){
  covars=1:max(unlist(ind.betas))
  if(is.null(nmodels)) nmodels=length(ind.betas)
  prob.covars=covars*0
  for(i in 1:nmodels) prob.covars[ind.betas[[i]]]=prob.covars[ind.betas[[i]]]+probs[i+1]
  if(any(prob.covars>0.5)){
    tt=which(prob.covars>0.5)
    res=which(laply(ind.betas,function(x) all(x%in%tt)&all(tt%in%x)))
  }else{res=0}
  return(res+1)
}


######## Main Function ########
# Calculate the BF, specified in typebf for all models.
# Return:
# - log of BF of 1 against the null;
# - model posterior probability;
# - HPPM model;
# - MPPM model;
# - Mean model size
# - Posterior Expected Model Size
allbf=function(ind.betas,dat,Bs,Bmode,parms.init,nsimul,use.pericchi2005.fbf,
               nprocess,log.post,typebf,force.L=NULL,use.MCMC=TRUE,log.bf.all=NULL,true.sim=NULL,
               use.true.sim=FALSE,logbf=NULL,no.bergscott.prior){
   nmodels=length(ind.betas)
  if(is.null(logbf)){
    tmp=mclapply(as.list(1:nmodels),function(ii){
      data0=as.matrix(cbind(dat[1:2],1))
      data1=as.matrix(cbind(dat[1:2],model.matrix(~.,data=dat[ind.betas[[ii]]+2])))
      typebf(data0,data1,parms.init.used=list(parms.init[[1]],parms.init[[ii+1]]),
             nsimul=nsimul,log.post=log.post,Bs=Bs[[ii]],Bmode=Bmode[[ii]],
             use.pericchi2005.fbf=use.pericchi2005.fbf,
             force.L=force.L,use.MCMC=use.MCMC,log.bf.all=log.bf.all[ii])},
                 mc.cores=nprocess)
    tmp=unlist(tmp)
    tmp=matrix(tmp,nmodels,length(tmp)/nmodels,byrow=TRUE,dimnames=list(NULL,names(tmp)[1:(length(tmp)/nmodels)]))
    logbf=NULL
    for(i in 1:ncol(tmp)) logbf=cbind(logbf,unlist(tmp[,i]))
    logbf=rbind(rep(0,ncol(tmp)),logbf)
    colnames(logbf)=colnames(tmp)
  }
  rownames(logbf)=c("0",paste(ind.betas))
  model.post=function(x,no.bergscott.prior){
    if(no.bergscott.prior){
      res=exp(x)/sum(exp(x))}else{
        ncov=c(0,laply(ind.betas,length))
	uncov=unique(ncov)
        p=max(ncov)
	redprior=1/(p+1)
	redprior=redprior/c(table(ncov))
	modprior=redprior[ncov+1]
        res=(exp(x)*modprior)/sum(exp(x)*modprior)
      }
    return(res)
  }
  ppost=cbind(aaply(logbf,2,model.post,no.bergscott.prior=no.bergscott.prior))
  if(nrow(ppost)<ncol(ppost)) ppost=t(ppost)
  if(ncol(ppost)>1){
    hppm.num=aaply(ppost,2,function(x) which(x==max(x))[1])}else{
      hppm.num=which(ppost==max(ppost))}
  if(ncol(ppost)>1){mppm.num=aaply(ppost,2,mppm.func,ind.betas=ind.betas,nmodels=nmodels)}else{
    mppm.num=mppm.func(ppost,ind.betas,nmodels)}
  hppm=aaply(hppm.num,1,function(x) rownames(logbf)[x])
  mppm=aaply(mppm.num,1,function(x) rownames(logbf)[x])
  pems=c(0,laply(ind.betas,function(x) length(x)))
  if(ncol(ppost)>1){pems=t(ppost)%*%pems}else{pems=sum(ppost*pems)}
  name.cov=colnames(dat)[-(1:2)]
  rownames(logbf)=c("*NULL*",aaply(rownames(logbf)[-1],1,function(x) paste(name.cov[eval(parse(text=x))],collapse="-")))
  rownames(ppost)=c("*NULL*",aaply(rownames(ppost)[-1],1,function(x) paste(name.cov[eval(parse(text=x))],collapse="-")))
  hppm=laply(hppm,function(x) paste(name.cov[eval(parse(text=x))],collapse="-"))
  hppm[hppm==""]="*NULL*"  
  mppm=laply(mppm,function(x) paste(name.cov[eval(parse(text=x))],collapse="-"))
  mppm[mppm==""]="*NULL*"
  if(length(hppm)>1) names(hppm)=names(mppm)=colnames(tmp)
  if(is.null(true.sim)){
    res.allbf=list(logbf=logbf,ppost=ppost,hppm=hppm,mppm=mppm,pems=pems)}else{
      ppost=as.matrix(ppost)
      ppost.true=ppost[true.sim,]
      if(length(ppost.true)==1) names(ppost.true)=NULL
      pems=as.matrix(pems)
      names(pems)=names(ppost.true)
      hppm.succ=hppm.num==true.sim
      mppm.succ=mppm.num==true.sim
      names(mppm.succ)=NULL
      res.allbf=list(ppost=ppost.true,hppm=hppm.succ,mppm=mppm.succ,pems=pems)
    }
  return(res.allbf)
}

analyse.data=function(datos,calc.cpbf=TRUE,calc.bic=TRUE,calc.fbf=TRUE,calc.fbfmode=TRUE,use.pericchi2005.fbf=TRUE,
                       calc.ibf=TRUE,nsimul=1000,force.L=TRUE,useB="Bmode",useB.ibf="Bmode",
                       use.Bsuni=FALSE,use.MCMC=FALSE,model,nprocess,true.sim=NULL,use.true.sim=FALSE,onlyone=FALSE,
                       no.bergscott.prior=TRUE,calc.ibfkm=FALSE){
  if(lmodel=="log.post.gamma") log.post=log.post.ggamma
  if(lmodel=="log.post.lnormal") log.post=log.post.lnormal
  if(lmodel=="log.post.weibull") log.post=log.post.weibull
  if(lmodel=="log.post.lnormal.cp") log.post=log.post.lnormal.cp
  
  # Model Index
  numcov=ncol(datos)-2
  ind.betas=index.models(numcov,onlyone=onlyone)
  n=nrow(datos)
  #if(use.true.sim){
    beta.trues=llply(ind.betas,function(x){
      rr=rep(0,numcov)
      rr[x]=1
      return(rr)})
#}
  # Parms init
  data0=as.matrix(cbind(datos[1:2],1))
  if(use.true.sim){tt1=c(0,0)}else{
    tt1=optim(par=c(0,0),fn=log.post,control=list(fnscale=-1),Data=data0,b=c(1,1),method=optimmethod)$par
  }
  tmpmyf=function(ii){
    data1=as.matrix(cbind(datos[1:2],model.matrix(~.,data=datos[ind.betas[[ii]]+2])))
    if(use.true.sim){pars=c(tt1,unlist(beta.trues[ii])[ind.betas[[ii]]])
                     logbf=log.post(parms=pars,Data=data1,b=c(1,1))-log.post(parms=tt1,Data=data0,b=c(1,1))
    }else{
       
      nunc=nunc.func(data1[,2])
      XtX.inv=func.XtX.inv(rel=data1[,2],XM=cbind(data1[,-(1:3)]),y=data1[,1])
      pars=optim(par=c(rep(0,2),unlist(beta.trues[ii])[ind.betas[[ii]]]),fn=log.post,control=list(fnscale=-1),Data=data1,b=c(1,1),nunc=nunc,XtX.inv=XtX.inv,method=optimmethod)$par
      logbf=log.bf.parziale(data0,data1,nsimul,parms.init.used=list(tt1,pars),b=c(1,1),log.post=log.post,use.MCMC=use.MCMC)
         
    }
    return(list(pars=pars,logbf=logbf))}
  ltt=alply(1:length(ind.betas),1,tmpmyf,.parallel=TRUE)
  parms.init=c(list(tt1))
  log.bf.all=NULL
  
    for(i in 1:length(ind.betas)){
      parms.init=c(parms.init,ltt[[i]][1])
      log.bf.all=c(log.bf.all,ltt[[i]][2])  
    }
 
  log.bf.all=unlist(log.bf.all)
  
  nummaxpar=max(unlist(llply(ind.betas,length)))+2
  # Calculates of Fraction B
  if(use.Bsuni){
    # All same B and same L
    allB=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=nummaxpar,use="all"))
    Bmode=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=nummaxpar,use="mode"))
    #Bmedian=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=numcov+2,use="median"))
  }else{
    # Different B and Different L
    allB=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=length(x)+2,use="all"))
    Bmode=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=length(x)+2,use="mode"))
    #Bmedian=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=length(x)+2,use="median"))
  }
  Bs=get(useB)
  
  #BF with Conventional Prior:
  if(calc.cpbf){
    res.cpbf=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                   nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                   typebf=NULL,use.MCMC=NULL,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior,logbf=t(t(c(0,log.bf.all))))
  }else{res.cpbf=NULL}
  
  # Fractional BF
  if(calc.fbf){
    res.fbf=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=use.pericchi2005.fbf,
                  typebf=fbf.ind.beta,use.MCMC=use.MCMC,log.bf.all=log.bf.all,true.sim=true.sim,
                  use.true.sim=use.true.sim,no.bergscott.prior=no.bergscott.prior)
    
  }else{res.fbf=NULL}
  # Intrinsic BF
  if(calc.ibf){
    num.max.fractions=NULL
    if(force.L) num.max.fractions=max(aaply(1:length(ind.betas),.margins=1,.fun=function(i) length(allB[[i]]$B[,1])))
    Bs=get(useB.ibf)
    res.ibf=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                  typebf=ibf.ind.beta,force.L=num.max.fractions,use.MCMC=use.MCMC,
                  log.bf.all=log.bf.all,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior)
  }else{res.ibf=NULL}
  # KM Intrinsic BF
  if(calc.ibfkm){
    num.max.fractions=NULL
    if(force.L) num.max.fractions=max(aaply(1:length(ind.betas),.margins=1,.fun=function(i) length(allB[[i]]$B[,1])))
    Bs=get(useB.ibf)
    res.ibfkm=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                    nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                    typebf=ibfkm.ind.beta,force.L=num.max.fractions,use.MCMC=use.MCMC,
                    log.bf.all=log.bf.all,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior)
  }else{res.ibfkm=NULL}
  # BIC
  if(calc.bic){
    res.bic=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                  typebf=bic.ind.beta,use.MCMC=NULL,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior)
  }else{res.bic=NULL}
  
  return(list(fbf=res.fbf,ibf=res.ibf,bic=res.bic,ibfkm=res.ibfkm,cpbf=res.cpbf))
}




analyse.data1cov=function(datos,calc.cpbf=TRUE,calc.bic=TRUE,calc.fbf=TRUE,calc.fbfmode=TRUE,use.pericchi2005.fbf=TRUE,
                      calc.ibf=TRUE,nsimul=1000,force.L=TRUE,useB="Bmode",useB.ibf="Bmode",
                      use.Bsuni=FALSE,use.MCMC=FALSE,model,nprocess,true.sim=NULL,use.true.sim=FALSE,onlyone=FALSE,
                      no.bergscott.prior=TRUE,calc.ibfkm=FALSE){
  if(lmodel=="log.post.gamma") log.post=log.post.ggamma
  if(lmodel=="log.post.lnormal") log.post=log.post.lnormal
  if(lmodel=="log.post.weibull") log.post=log.post.weibull
  if(lmodel=="log.post.lnormal.cp") log.post=log.post.lnormal.cp
  
  # Model Index
  numcov=ncol(datos)-2
  ind.betas=index.models(numcov,onlyone=onlyone)
  n=nrow(datos)
  if(use.true.sim){
    beta.trues=llply(ind.betas,function(x){
      rr=rep(0,numcov)
      rr[x]=1
      return(rr)})}
  # Parms init
  data0=as.matrix(cbind(datos[1:2],1))
  if(use.true.sim){tt1=c(0,0)}else{
    tt1=optim(par=c(0,0),fn=log.post,control=list(fnscale=-1),Data=data0,b=c(1,1),method=optimmethod)$par
  }
  tmpmyf=function(ii){
    data1=as.matrix(cbind(datos[1:2],model.matrix(~.,data=datos[ind.betas[[ii]]+2])))
    if(use.true.sim){pars=c(tt1,unlist(beta.trues[ii])[ind.betas[[ii]]])
                     logbf=log.post(parms=pars,Data=data1,b=c(1,1))-log.post(parms=tt1,Data=data0,b=c(1,1))
    }else{
      nunc=nunc.func(data1[,2])
      XtX.inv=func.XtX.inv(rel=data1[,2],XM=cbind(data1[,-(1:3)]),y=data1[,2])
      pars=optim(par=c(0,rep(0,ncol(data1)-2)),fn=log.post,control=list(fnscale=-1),Data=data1,b=c(1,1),nunc=nunc,XtX.inv=XtX.inv,method=optimmethod)$par
      logbf=log.bf.parziale(data0,data1,nsimul,parms.init.used=list(tt1,pars),b=c(1,1),log.post=log.post,use.MCMC=use.MCMC)
    }
    return(list(pars=pars,logbf=logbf))}
  ltt=tmpmyf(1)
  parms.init=c(list(tt1))
  log.bf.all=NULL
  
  parms.init=ltt[1]
  log.bf.all=ltt[2] 
  
   log.bf.all=unlist(log.bf.all)
  
  nummaxpar=max(unlist(llply(ind.betas,length)))+2
  # Calculates of Fraction B
  if(use.Bsuni){
    # All same B and same L
    allB=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=nummaxpar,use="all"))
    Bmode=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=nummaxpar,use="mode"))
    #Bmedian=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=numcov+2,use="median"))
  }else{
    # Different B and Different L
    allB=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=length(x)+2,use="all"))
    Bmode=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=length(x)+2,use="mode"))
    #Bmedian=llply(ind.betas,function(x) prob.b(n=n,ncens=n-sum(datos[,2]),p=length(x)+2,use="median"))
  }
  Bs=get(useB)
  
  #BF con Conventional Prior:
  if(calc.cpbf){
    res.cpbf=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                   nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                   typebf=NULL,use.MCMC=NULL,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior,logbf=t(t(c(0,log.bf.all))))
  }else{res.cpbf=NULL}
  
  # Fractional BF
  if(calc.fbf){
    res.fbf=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=use.pericchi2005.fbf,
                  typebf=fbf.ind.beta,use.MCMC=use.MCMC,log.bf.all=log.bf.all,true.sim=true.sim,
                  use.true.sim=use.true.sim,no.bergscott.prior=no.bergscott.prior)
    
  }else{res.fbf=NULL}
  # Intrinsic BF
  if(calc.ibf){
    num.max.fractions=NULL
    if(force.L) num.max.fractions=max(aaply(1:length(ind.betas),.margins=1,.fun=function(i) length(allB[[i]]$B[,1])))
    Bs=get(useB.ibf)
    res.ibf=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                  typebf=ibf.ind.beta,force.L=num.max.fractions,use.MCMC=use.MCMC,
                  log.bf.all=log.bf.all,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior)
  }else{res.ibf=NULL}
  # KM Intrinsic BF
  if(calc.ibfkm){
    num.max.fractions=NULL
    if(force.L) num.max.fractions=max(aaply(1:length(ind.betas),.margins=1,.fun=function(i) length(allB[[i]]$B[,1])))
    Bs=get(useB.ibf)
    res.ibfkm=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                  typebf=ibfkm.ind.beta,force.L=num.max.fractions,use.MCMC=use.MCMC,
                  log.bf.all=log.bf.all,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior)
  }else{res.ibfkm=NULL}
  # BIC
  if(calc.bic){
    res.bic=allbf(ind.betas=ind.betas,dat=datos,Bs=Bs,Bmode=Bmode,parms.init=parms.init,
                  nsimul=nsimul,nprocess=nprocess,log.post=log.post,use.pericchi2005.fbf=NULL,
                  typebf=bic.ind.beta,use.MCMC=NULL,true.sim=true.sim,no.bergscott.prior=no.bergscott.prior)
  }else{res.bic=NULL}
 
  return(list(fbf=res.fbf,ibf=res.ibf,bic=res.bic,ibfkm=res.ibfkm,cpbf=res.cpbf))
}


### Only BIC and FBF from pericchi 2005 approximation
analyse.data.bic.fbf=function(datos,Bs,Bmode,model,nprocess,true.sim=NULL){
  if(lmodel=="log.post.lnormal") log.post=log.post.lnormal
  if(lmodel=="log.post.weibull") log.post=log.post.weibull
  
  # Calculus of Maximum(S) of Log Likelihoods for each model and each possible fraction
  # For the null model
  data0=as.matrix(cbind(datos[1:2],1))
  log0f=rep(NA,nrow(allB[[1]]$B))
  for(i in 1:nrow(allB[[1]]$B)) log0f[i]=log.marg.laplace(log.post, parms.hat=c(0,0), data0, b=allB[[1]]$B[i,])
  log0.N=rep(log.marg.laplace(log.post, parms.hat=c(0,0), data0, b=c(1,1)),nrow(allB[[1]]$B))

  # For the rest of models
  tts=mcmapply(function(ii){
    data1=as.matrix(cbind(datos[1:2],model.matrix(~.,data=datos[ind.betas[[ii]]+2])))
    inits=c(0,rep(0,ncol(data1)-2))
    nfrac=nrow(allB[[ii+1]]$B)
    logj.N=rep(log.marg.laplace(log.post,parms.hat=inits,dat=data1,b=c(1,1)),nfrac)
    logjf=logj0f=rep(NA,nfrac)
    for(i in 1:nfrac){
      logjf[i]=log.marg.laplace(log.post,parms.hat=inits,dat=data1,b=allB[[ii+1]]$B[i,])
      logj0f[i]=log.marg.laplace(log.post,parms.hat=c(0,0),dat=data0,b=allB[[ii+1]]$B[i,])
    }
      return(list(logjf=logjf,logj0f=logj0f,logj.N=logj.N))},
               1:length(ind.betas),mc.cores=nprocess)
  
  llsf=c(list(log0f),tts[1,])
  lls0f=c(list(log0f),tts[2,])
  lls.N=c(list(log0.N),tts[3,])
  
  # Calculates of Maximum of Log Likelihoods and Model dimension
  data0=as.matrix(cbind(datos[1:2],1))
  log0=optim(par=c(0,0),fn=log.post,control=list(fnscale=-1),Data=data0,b=c(1,1),method=optimmethod)$value
  np0=2
  tts=mcmapply(function(ii){
    data1=as.matrix(cbind(datos[1:2],model.matrix(~.,data=datos[ind.betas[[ii]]+2])))
    inits=c(0,rep(0,ncol(data1)-2))
    logj=optim(par=inits,fn=log.post,control=list(fnscale=-1),Data=data1,b=c(1,1))$value
    return(list(logj=logj,npj=length(inits)))},1:length(ind.betas),mc.cores=nprocess)
  lls=c(log0,unlist(tts[1,]))
  npj=c(np0,unlist(tts[2,]))
  
  #logRf=log0f;npRf=npjf[[1]] # Encompassing from BELOW
  logR=log0;npR=np0 # Encompassing from BELOW
  # logR=lls[length(lls)];npR=npj[length(npj)] # Encompassing from ABOVE
  
  # BIC
  # Likelihood Evidence 
  EL=(lls-logR)
  # Difference of Dimension
  dd=(npj-npR)
  # Log of BIC
  lbic=EL-log(nunc)*dd*0.5
  
  nmodels=length(llsf)
  lmfbf=lfbfmode=rep(NA,nmodels)
  for(i in 1:length(llsf)){
#    allfbf=(lls.N[[i]]-lls.N[[nmodels]])-(llsf[[i]]-llsf[[nmodels]])
    allfbf=(lls.N[[i]]-lls.N[[1]])-(llsf[[i]]-lls0f[[i]])
    lmfbf[i]=log(sum(exp(allfbf)*Bs[[i]]$B.prob))
#    lmfbf[i]=sum(allfbf*Bs[[i]]$B.prob)
    lfbfmode[i]=allfbf[Bs[[i]]$B.prob==max(Bs[[i]]$B.prob)]
  }

  logbf=cbind(lbic,lfbfmode,lmfbf)
  
  colnames(logbf)=c("bic","fbfmode","mfbf")
  res=allbf(ind.betas,dat=datos,Bs=NULL,Bmode=NULL,parms.init=NULL,nsimul=NULL,use.pericchi2005.fbf=NULL,
            nprocess=NULL,log.post=NULL,typebf=NULL,force.L=NULL,use.MCMC=NULL,log.bf.all=NULL,
            true.sim=true.sim,use.true.sim=NULL,logbf=logbf)
  return(res)
}



############################ ANALYSE SIMULATIONS ######################
collect.sim=function(results,ibfkm=FALSE){
  aibflab=c("aibf","kmaibf")[1+ibfkm]
  mibflab=c("mibf","kmmibf")[1+ibfkm]
  write.table(x=results,file="tt.csv",sep=";",quote=FALSE)
  res=read.table("tt.csv",sep=";")
  res=data.frame(res[,-(1:2)])
  rownames(res)=NULL
  nsim=nrow(res)
  namcol=colnames(res)
  ids=grep("ppost",namcol)
  repsim=rep(1:nsim,length(ids))
  ppost=unlist(res[ids])
  typebf=namcol[ids]
  typebf[grep("mfbf",typebf)]="mfbf"
  typebf[grep("cpbf",typebf)]="cpbf"
  typebf[grep("fbfmode",typebf)]="fbfmo"
  typebf[grep("bic",typebf)]="bic"
  typebf[grep("aibf",typebf)]=aibflab
  typebf[grep("mibf",typebf)]=mibflab
  typebf=factor(typebf[gl(length(ids),nsim)],levels=c(aibflab,mibflab,"mfbf","fbfmo","bic","cpbf"))
  hppm=unlist(res[grep("hppm",namcol)])
  mppm=unlist(res[grep("mppm",namcol)])
  allsim=res[repsim,1:11]
  allsim=cbind(allsim,typebf,hppm,mppm,ppost)
  return(allsim)
}

calc.pems=function(results){
  pems=data.frame(results[grep("pems",colnames(results))],stringsAsFactors=FALSE)
  n=as.numeric(as.character(results$sim.id.n,pems))
  for(i in 1:ncol(pems)) pems[,i]=as.numeric(as.character(pems[,i]))
  tt=as.numeric(as.character(unique(results$sim.id.full.model.size)))
  if(tt>15){onlyone=TRUE}else{onlyone=FALSE}
  model.ids=c(0,index.models(tt,onlyone=onlyone))
  model.trues=laply(model.ids,function(x) length(x)*(min(x)>0))
  msize=model.trues[as.numeric(as.character(results$sim.id.true.model))]
  colnames(pems)[grep("mfbf",colnames(pems))]="mfbf"
  colnames(pems)[grep("cpbf",colnames(pems))]="cpbf"
  colnames(pems)[grep("fbfmode",colnames(pems))]="fbfmo"
  colnames(pems)[grep("bic",colnames(pems))]="bic"
  colnames(pems)[grep("aibf",colnames(pems))]="aibf"
  colnames(pems)[grep("mibf",colnames(pems))]="mibf"
  typebf=colnames(pems)[gl(ncol(pems),nrow(pems))]
  typebf=factor(typebf,levels=c("aibf","mibf","mfbf","fbfmo","bic","cpbf"))
  msize=rep(msize,ncol(pems))
  n=rep(n,ncol(pems))
  n=factor(n,levels=sort(unique(n)))
  tt=NULL
  for(i in 1:ncol(pems)) tt=c(tt,pems[,i])
  pems=tt
  return(data.frame(pems,n,msize,typebf))
}


analyse.res.dataset=function(fbf,bic,ibf,nmodel=2,decnn=3){
  ppmfbf=round(fbf$ppost[order(fbf$ppost[,1],decreasing=TRUE)[1:nmodel],1],decnn)
  ppfbfmode=round(fbf$ppost[order(fbf$ppost[,2],decreasing=TRUE)[1:nmodel],2],decnn)
  ppbic=round(bic$ppost[order(bic$ppost,decreasing=TRUE)[1:nmodel],],decnn)
  ppaibf=round(ibf$ppost[order(ibf$ppost[,1],decreasing=TRUE)[1:nmodel],1],decnn)
  ppmibf=round(ibf$ppost[order(ibf$ppost[,2],decreasing=TRUE)[1:nmodel],2],decnn)
  allprobs=sort(c(ppmfbf,ppfbfmode,ppbic,ppaibf,ppmibf),decreasing=TRUE)
  allmodels=unique(names(allprobs))
  ppdata=array(NA,dim=c(length(allmodels),7))
  ppdata=data.frame(ppdata)
  colnames(ppdata)=c("mFBF","FBFmo","BIC","aibf","mibf","HPPM","MPPM")
  rownames(ppdata)=allmodels
  allhppm=c(fbf$hppm[1],fbf$hppm[2],bic$hppm,ibf$hppm[1],ibf$hppm[2])
  allmppm=c(fbf$mppm[1],fbf$mppm[2],bic$mppm,ibf$mppm[1],ibf$mppm[2])
  i=1
  for(mm in allmodels){
    tmphppm=paste(colnames(ppdata)[which(allhppm==mm)],collapse=",")
    tmpmppm=paste(colnames(ppdata)[which(allmppm==mm)],collapse=",")
    ttprobs=c(ifelse(length(ppmfbf[names(ppmfbf)==mm])>0,ppmfbf[names(ppmfbf)==mm],NA),
              ifelse(length(ppfbfmode[names(ppfbfmode)==mm])>0,ppfbfmode[names(ppfbfmode)==mm],NA),
              ifelse(length(ppbic[names(ppbic)==mm])>0,ppbic[names(ppbic)==mm],NA),
              ifelse(length(ppaibf[names(ppaibf)==mm])>0,ppaibf[names(ppaibf)==mm],NA),
              ifelse(length(ppmibf[names(ppmibf)==mm])>0,ppmibf[names(ppmibf)==mm],NA),
              ifelse(length(tmphppm)>0,tmphppm,NA),ifelse(length(tmpmppm)>0,tmpmppm,NA))
    ppdata[i,]=ttprobs
    i=i+1
  }
  pems=c(round(c(fbf$pems,bic$pems,ibf$pems),1),"","")
  ppdata=rbind(ppdata,pems)
  rownames(ppdata)[nrow(ppdata)]="PEMS"
  ppdata[is.na(ppdata)]=""
  #xtable(ppdata,align="l|ccccccc",caption="Results for XXX data",label="tab-xxx")
  return(ppdata)
}

aggiusta.results=function(results){
  ppost=rep(NA,nrow(results))
  true.mod.ids=as.numeric(as.character(results$sim.id.true.model))
  for(i in 1:nrow(results)) ppost[i]=as.numeric(as.character(results[[21+true.mod.ids[i]]][i]))
  results=results[-(22:29)]
  results=data.frame(results,cpbf.ppost=ppost)
  
  true.mod.lab=c("*NULL*","x1","x2,","x3,","x1-x2","x1-x3,","x2-x3,","x1-x2-x3")[true.mod.ids]
  results$cpbf.hppm.1=results$cpbf.hppm.1==true.mod.lab
  results$cpbf.mppm.1=results$cpbf.mppm.1==true.mod.lab
  return(results)
}


