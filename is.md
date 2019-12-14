Intro
=====

This code implements, the importance sampling approach to obtain
marginals and hence Bayes Factors.

This does not implement stochastic model space exploration and all
possible models (2<sup>*k*</sup>, where *k* is the number of covariates)
are visited.

Data
====

We consider the heart transplant survival dataset and two spurious
uncorrelated variable. This is just one draw for the simulation study
described in the paper.

    rm(list=ls()) 
    set.seed(17) 
    library(compiler) 
    library(mvtnorm) 
    library(plyr) 
    library(doParallel) 
    source("library-bf-cens-conv-prior.R") 
    source("library-bf-cens-conv-prior-2.R") 
    library(survival) 
    data(heart) 
    jasa=jasa[jasa$transplant==1,] 
    surv.time=jasa$fu.date-jasa$accept.dt
    fecha.fin=as.Date(0, origin = "1974-04-01") 
    cens.time=fecha.fin-jasa$accept.dt  
    cens.time[jasa$fustat==0]=surv.time[jasa$fustat==0] 
    dat.heart=data.frame(as.numeric(surv.time),as.numeric(jasa$fustat),as.numeric(cens.time),jasa$age,rnorm(nrow(jasa)),rnorm(nrow(jasa))) 
    colnames(dat.heart)=c("survival","rel","cens","age","spur1","spur2") 
    dat.heart=dat.heart[dat.heart$survival>0,] 
    summary(dat.heart)

    ##     survival           rel              cens           age       
    ##  Min.   :   4.0   Min.   :0.0000   Min.   :  38   Min.   :19.55  
    ##  1st Qu.:  71.0   1st Qu.:0.0000   1st Qu.: 487   1st Qu.:42.50  
    ##  Median : 206.0   Median :1.0000   Median : 941   Median :47.99  
    ##  Mean   : 414.5   Mean   :0.6522   Mean   :1060   Mean   :46.03  
    ##  3rd Qu.: 619.0   3rd Qu.:1.0000   3rd Qu.:1676   3rd Qu.:52.03  
    ##  Max.   :1799.0   Max.   :1.0000   Max.   :2277   Max.   :64.41  
    ##      spur1             spur2         
    ##  Min.   :-3.5561   Min.   :-3.32032  
    ##  1st Qu.:-0.6166   1st Qu.:-0.64045  
    ##  Median : 0.1510   Median :-0.04963  
    ##  Mean   : 0.1375   Mean   :-0.14477  
    ##  3rd Qu.: 0.7721   3rd Qu.: 0.59098  
    ##  Max.   : 2.4423   Max.   : 2.10842

Analysis
========

The null model
--------------

Data for the null model are

    data0=as.matrix(cbind(log(dat.heart$survival),dat.heart$rel,1))

IS approximation for the marginal

    tt0=optim(par=c(0,0),fn=log.post.lnormal,control=list(fnscale=-1),Data=data0,method=optimmethod,hessian=TRUE)
    var.prop=-solve(tt0$hessian)

    nsimul=10000
    parms.hat0=tt0$par
    res=post(log.post.lnormal,nsimul=nsimul,parms.init=parms.hat0,dat=data0)
    res$par=res$par[-(1:(nsimul*0.3)),]
    mu.prop=apply(res$par,2,mean)
    Sigma.prop=cov(res$par)
    nsim=40000
    thetas=rmvt(nsim, delta = mu.prop, sigma = Sigma.prop,df=3)

    lnum=apply(thetas,1,log.post.lnormal.0.our,Data=data0)
    lden=apply(thetas,1,dmvt,delta=mu.prop, sigma =Sigma.prop,df=3,log=TRUE)

    (lm0=log(mean(exp(lnum-lden))))

    ## [1] -71.06686

Rest of the models
------------------

The rest of the models is defined upon the full matrix. All models will
be examined (no random search).

    var.name=c("age","spur1","spur2") 
    k=length(var.name)
    n=nrow(dat.heart)
    Xfull=dat.heart[,var.name]
    Xfull=scale(Xfull) 
    colnames(Xfull)=var.name 
    Xfull=as.data.frame(Xfull) 
    mod.list=index.models(ncol(Xfull)) 
    nmodels=length(mod.list) 

Parallel calculation:

    cl <- makePSOCKcluster(4) 
    registerDoParallel(cl) 

These are the examined models:

    llply(mod.list,function(x) colnames(Xfull)[x]) 

    ## [[1]]
    ## [1] "age"
    ## 
    ## [[2]]
    ## [1] "spur1"
    ## 
    ## [[3]]
    ## [1] "spur2"
    ## 
    ## [[4]]
    ## [1] "age"   "spur1"
    ## 
    ## [[5]]
    ## [1] "age"   "spur2"
    ## 
    ## [[6]]
    ## [1] "spur1" "spur2"
    ## 
    ## [[7]]
    ## [1] "age"   "spur1" "spur2"

We have 7 models more than the null, whose marginals are approximated
with IS with *g* is integrated over the robust prior (uncomment ) or
inverse Gamma (uncomment ) or fixed (uncomment ) with this code.

    lm1 <- foreach(i=1:nmodels, .combine=c) %dopar% { 
        source("library-bf-cens-conv-prior.R") 
        source("library-bf-cens-conv-prior-2.R")
        library(mvtnorm)
        #robust# marg1.allpara.robust(y=log(dat.heart$survival),rel=dat.heart$rel,ct=log(dat.heart$cens),X=Xfull[mod.list[[i]]],k.star=0,nsim.marg = 10000,nsim.post = 5000)$marg
     #ig#
      marg1.allpara.igamma(y=log(dat.heart$survival),rel=dat.heart$rel,cens=log(dat.heart$cens),X=Xfull[mod.list[[i]]],k.star=0,nsim.marg = 10000,nsim.post = 5000)$marg 
      #gfixed#  marg1.allpara.gfixed(y=log(dat.heart$survival),rel=dat.heart$rel,cens=log(dat.heart$cens),X=Xfull[mod.list[[i]]],k.star=0,nsim.marg = 10000,nsim.post = 5000)$marg
    }
    stopCluster(cl)
    lm1=log(lm1)

BF, model posterior probabilities and posterior inclusion probabilities
=======================================================================

BFs and model posterior probabilities
-------------------------------------

    mdim=unlist(lapply(mod.list,function(x) length(x))) 
    norm.mdim=table(mdim) 
    mod.prior=1/norm.mdim[mdim] 
    mod.prior=c(mod.prior,1) 
    mod.prior=mod.prior/sum(mod.prior) 
    BF=exp(lm1-lm0) 
    mod.pp=prob.post(BF,mod.list,scott=TRUE) 

A matrix with the summary:

    results<- matrix(0, ncol=ncol(Xfull)+4, nrow=nmodels+1) 
    colnames(results)<- c(colnames(Xfull), "marg", "BF", "PriorProb", "PosteriorProb") 
    for (i in 1:nmodels){ 
      results[i,mod.list[[i]]]<- 1 
      results[i,ncol(Xfull)+(1:2)]<- c(lm1[i],BF[i]) 
      results[i,ncol(Xfull)+3]<- mod.prior[i] 
    }
    results[nmodels+1,ncol(Xfull)+(1:3)]<- c(lm0, 1,mod.prior[nmodels+1]) 
    results[,"PosteriorProb"]<- c(mod.pp[-1],mod.pp[1]) 
    rownames(results)=c(paste("Model ",1:nmodels),"Null Model") 
    results 

    ##            age spur1 spur2      marg         BF  PriorProb PosteriorProb
    ## Model  1     1     0     0 -71.45198 0.68037435 0.08333333    0.15814293
    ## Model  2     0     1     0 -73.10694 0.13001932 0.08333333    0.03022106
    ## Model  3     0     0     1 -73.28758 0.10853115 0.08333333    0.02522646
    ## Model  4     1     1     0 -73.06747 0.13525362 0.08333333    0.03143770
    ## Model  5     1     0     1 -73.12746 0.12737802 0.08333333    0.02960713
    ## Model  6     0     1     1 -74.86450 0.02242376 0.08333333    0.00521207
    ## Model  7     1     1     1 -74.48526 0.03276486 0.25000000    0.02284712
    ## Null Model   0     0     0 -71.06686 1.00000000 0.25000000    0.69730553

High posterior model
--------------------

    hpm=order(results[,"PosteriorProb"],decreasing = TRUE)[1] 
    rownames(results)[hpm] 

    ## [1] "Null Model"

Posterior inclusion probabilities and median probability model
--------------------------------------------------------------

    (pip=colSums(t(t(results[,1:ncol(Xfull)]*results[,"PosteriorProb"])))) 

    ##        age      spur1      spur2 
    ## 0.24203488 0.08971795 0.08289278

    mpm=(pip>0.5)*1 
    mpm=apply(results[,1:k],1,function(x) all(x==mpm)) 
    rownames(results)[mpm] 

    ## [1] "Null Model"
