A Model selection approach for Variable selection with censored data (R
code)
================
María Eugenia Castellanos, Gonzalo García Donato and Stefano Cabras

# Intro

This code implements, for each entertained model
\(\Sigma^M_\gamma(\hat\beta_0,\hat\sigma)\) as the prior covariance
matrix, where the point estimators are the posterior means of the
corresponding parameters under the null model.

Moreover, all marginals are approximated with Laplace.

This does not implement stochastic model space exploration and all
possible models (\(2^k\), where \(k\) is the number of covariates) are
visited.

# Data

We consider the heart transplant survival dataset and two spuroius
uncorrelated variable. This is just one draw for the simulation study
described in the paper.

``` r
rm(list=ls())
library(compiler)
library(mvtnorm)
library(plyr)
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` r
source("library-bf-cens-conv-prior.R")
```

    ## Loading required package: coda

    ## Loading required package: MASS

    ## ##
    ## ## Markov Chain Monte Carlo Package (MCMCpack)

    ## ## Copyright (C) 2003-2019 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park

    ## ##
    ## ## Support provided by the U.S. National Science Foundation

    ## ## (Grants SES-0350646 and SES-0350613)
    ## ##

    ## 
    ## Attaching package: 'LearnBayes'

    ## The following object is masked from 'package:MCMCpack':
    ## 
    ##     rdirichlet

    ## dummies-1.5.6 provided by Decision Patterns

``` r
source("library-bf-cens-conv-prior-2.R")
library(survival)
data(heart)
#Seleccionamos solo los pacientes que han recibido trasplante
jasa=jasa[jasa$transplant==1,]
#El tiempo de seguimiento es:
surv.time=jasa$fu.date-jasa$accept.dt  
#La fecha a la que se finaliza el estudio es:
fecha.fin=as.Date(0, origin = "1974-04-01")
cens.time=fecha.fin-jasa$accept.dt 
#Para los censurados, nos aseguramos de definir bien el tiempo de supervivencia
cens.time[jasa$fustat==0]=surv.time[jasa$fustat==0]
#Usando edad:
dat.heart=data.frame(as.numeric(surv.time),as.numeric(jasa$fustat),as.numeric(cens.time),jasa$age,rnorm(nrow(jasa)),rnorm(nrow(jasa)))
colnames(dat.heart)=c("survival","rel","cens","age","spur1","spur2")

dat.heart=dat.heart[dat.heart$survival>0,]
summary(dat.heart)
```

    ##     survival           rel              cens           age       
    ##  Min.   :   4.0   Min.   :0.0000   Min.   :  38   Min.   :19.55  
    ##  1st Qu.:  71.0   1st Qu.:0.0000   1st Qu.: 487   1st Qu.:42.50  
    ##  Median : 206.0   Median :1.0000   Median : 941   Median :47.99  
    ##  Mean   : 414.5   Mean   :0.6522   Mean   :1060   Mean   :46.03  
    ##  3rd Qu.: 619.0   3rd Qu.:1.0000   3rd Qu.:1676   3rd Qu.:52.03  
    ##  Max.   :1799.0   Max.   :1.0000   Max.   :2277   Max.   :64.41  
    ##      spur1             spur2         
    ##  Min.   :-1.9056   Min.   :-2.31171  
    ##  1st Qu.:-0.7565   1st Qu.:-0.76695  
    ##  Median : 0.3227   Median :-0.03600  
    ##  Mean   : 0.1096   Mean   :-0.03167  
    ##  3rd Qu.: 0.9251   3rd Qu.: 0.57404  
    ##  Max.   : 1.5054   Max.   : 2.56103

# Analisis

## The null model

Data for the null model is

``` r
data0=as.matrix(cbind(log(dat.heart$survival),dat.heart$rel,1))
```

Laplace approximation for the marginal of the null model and calcualtion
of
\(\Sigma^M_\gamma(\hat\beta_0,\hat\sigma)\):

``` r
tt0=optim(par=c(0,0),fn=log.post.lnormal,control=list(fnscale=-1),Data=data0,method=optimmethod,hessian=TRUE)
var.prop=-solve(tt0$hessian)

parms.hat0=tt0$par
beta0.hat=parms.hat0[2]
lsigma.hat=parms.hat0[1]

hessi0=-tt0$hessian
map.coef=as.vector(tt0$par)
p=dim(hessi0)[1]

lm0=log.post.lnormal(parms=map.coef,Data=data0)+p/2*log(2*pi)-0.5*log(det(hessi0))
```

## Rest of the models

The rest of the models is defined upon the full matrix. All models will
be examined (no random search).

``` r
var.name=c("age","spur1","spur2")
k=length(var.name)
n=nrow(dat.heart)
Xfull=dat.heart[,var.name]
#Centramos la matriz X.full:
Xfull=scale(Xfull)
colnames(Xfull)=var.name
Xfull=as.data.frame(Xfull)
mod.list=index.models(ncol(Xfull))
nmodels=length(mod.list)
```

If parallel calculation is needed, please uncomment this:

``` r
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
```

These are the examined models:

``` r
llply(mod.list,function(x) colnames(Xfull)[x])
```

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
with Laplace with this code.

``` r
lm1 <- foreach(i=1:nmodels, .combine=c) %dopar% {
    source("library-bf-cens-conv-prior.R")
    source("library-bf-cens-conv-prior-2.R")
    marg1.allpara.gfixed.laplace(y=log(dat.heart$survival),
                                rel=dat.heart$rel,ct=log(dat.heart$cens),
                                X=as.matrix(Xfull[mod.list[[i]]]),
                                k.star=0,nsim.marg = NULL,
                                nsim.post = NULL,beta0.hat,lsigma.hat)$marg
}
```

# BF, model posterior probabilities and posterior inclusion probabilities

## BFs and model posterior probabilities

``` r
mdim=unlist(lapply(mod.list,function(x) length(x)))
norm.mdim=table(mdim)
mod.prior=1/norm.mdim[mdim]
mod.prior=c(mod.prior,1)
mod.prior=mod.prior/sum(mod.prior)
BF=exp(lm1-lm0)
mod.pp=prob.post(BF,mod.list,scott=TRUE)
#A matrix with the summary:
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
```

    ##            age spur1 spur2      marg         BF  PriorProb PosteriorProb
    ## Model  1     1     0     0 -71.22097 0.87889172 0.08333333   0.194261757
    ## Model  2     0     1     0 -72.87623 0.16790611 0.08333333   0.037112349
    ## Model  3     0     0     1 -73.07728 0.13732511 0.08333333   0.030353020
    ## Model  4     1     1     0 -73.04442 0.14191252 0.08333333   0.031366977
    ## Model  5     1     0     1 -73.21154 0.12007287 0.08333333   0.026539751
    ## Model  6     0     1     1 -74.90237 0.02213727 0.08333333   0.004893008
    ## Model  7     1     1     1 -75.07254 0.01867322 0.25000000   0.012382045
    ## Null Model   0     0     0 -71.09188 1.00000000 0.25000000   0.663091093

## High posterior model

``` r
hpm=order(results[,"PosteriorProb"],decreasing = TRUE)[1]
rownames(results)[hpm]
```

    ## [1] "Null Model"

## Posterior inclusion probabilities and median probability model

``` r
(pip=colSums(t(t(results[,1:ncol(Xfull)]*results[,"PosteriorProb"]))))
```

    ##        age      spur1      spur2 
    ## 0.26455053 0.08575438 0.07416782

``` r
mpm=(pip>0.5)*1
mpm=apply(results[,1:k],1,function(x) all(x==mpm))
rownames(results)[mpm]
```

    ## [1] "Null Model"
