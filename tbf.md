Intro
=====

This code implements, the analysis with test-based Bayes factors (TBF)
originally proposed by Johnson (2008) and Hu and Johnson (2009) later
revisited by Held, Sabanés Bové, and Gravestock (2015),

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

Posterior inclusion probabilities based on TBF
----------------------------------------------

The full matrix is

    var.name=c("age","spur1","spur2") 
    k=length(var.name)
    n=nrow(dat.heart)
    Xfull=dat.heart[,var.name]
    Xfull=scale(Xfull)
    colnames(Xfull)=var.name 
    Xfull=as.data.frame(Xfull) 
    mod.list=index.models(ncol(Xfull)) 
    nmodels=length(mod.list) 

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

These are the inclusion probabilities according to different versions of
the TBF:

    source("TBFfunctions.R")
    ipTBF(Xfull, y=dat.heart$survival, delta=dat.heart$rel, g.param="EB")

    ##       age     spur1     spur2 
    ## 0.5989054 0.4894859 0.4864764

    ipTBF(Xfull, y=dat.heart$survival, delta=dat.heart$rel, g.param="g=nu")

    ##       age     spur1     spur2 
    ## 0.3021753 0.1103147 0.1011396

    ipTBF(Xfull, y=dat.heart$survival, delta=dat.heart$rel, g.param="g=n")

    ##        age      spur1      spur2 
    ## 0.25252193 0.08415812 0.07625103

    ZSadapt.ipTBF(Xfull, y=dat.heart$survival, delta=dat.heart$rel) 

    ##       age     spur1     spur2 
    ## 0.2685626 0.1103536 0.1026478

References
==========

Held, L., D. Sabanés Bové, and I Gravestock. 2015. “Approximate Bayesian
Model Selection with the Deviance Statistic.” *Statistical Science* 30:
242–57.

Hu, J., and V.E. Johnson. 2009. “Bayesian Model Selection Using Test
Statistics.” *Journal of the Royal Statistical Society. Series B
(Methodological)* 71 (1): 143–58.

Johnson, V.E. 2008. “Properties of Bayes Factors Based on Test
Statistics.” *Scandinavian Journal of Statistics* 35 (2): 354–68.
