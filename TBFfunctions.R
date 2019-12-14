#Functions to compute inclusion probabilities within TBF's and hyper-TBF's

integer.base.b <-function(x, k){
        #k is the number of covariates in the most complex model
        ndigits <- (floor(logb(x, base=2))+1)
        res<- rep(0, ndigits)
        for(i in 1:ndigits){#i <- 1
                res[ndigits-i+1] <- (x %% 2)
                x <- (x %/% 2)
        }
        return(c(rep(0,k-ndigits),res))
}

#inclusion probs in the TBF
ipTBF<- function(Xfull, y, delta, g.param="EB"){
	k<- dim(Xfull)[2]; n<- dim(Xfull)[1]
	mBF<- matrix(0, nc=k+3, nr=2^k, byr=T); colnames(mBF)<- c(colnames(Xfull),"dev","tBF","postprob")
	#first the null:
	l0<- logLik(survreg(Surv(y, delta)~1, dist="lognormal"))
	mBF[2^k, k+1]<- 0; mBF[2^k, k+2]<- 1
	for (i in 1:(2^k-1)){
		#binary expression of the model
		mBF[i, 1:k]<- integer.base.b(i, k)
		#obtain the deviance
		fit<- survreg(Surv(y, delta)~., data=as.data.frame(cbind(y, delta, Xfull[, mBF[i, 1:k]==1])), dist="lognormal")
		#deviance:
		mBF[i, k+1]<- 2*(logLik(fit)-l0)		
		dj<- sum(mBF[i, 1:k])
	
		v<- 0
		if (g.param=="EB"){
			g<- max(mBF[i, k+1]/dj-1, 0); v<- 1
		}
		if (g.param=="g=n"){
		  g<- n; v<- 1
		}
		if (g.param=="g=nu"){
			g<-sum(delta); v<- 1
		}
		if (v==0) stop("Please specify a valid g.param\n")
	
		#test-based BF
		mBF[i, k+2]<- exp(-.5*dj*log(g+1)+g*mBF[i, k+1]/(2*(g+1)))		
	}
	  #posterior probs
		mBF[,k+3]<- mBF[, k+2]/((k+1)*choose(k, apply(mBF[,1:k], 1, sum)))
		mBF[,k+3]<- mBF[, k+3]/sum(mBF[,k+3])
	

	#inclusion probs:
	#(with fixed g=n)
	ip<- colSums(mBF[,1:k]*mBF[,k+3])
	#[1] 0.9843966 0.8179956 0.5906798 0.6298049
	#(with EB-basd g)
	return(colSums(mBF[,1:k]*mBF[,k+3]))
}


#inclusion probs in the hyper-g (adapted ZS) TBF
ZSadapt.ipTBF<- function(Xfull, y, delta){
	k<- dim(Xfull)[2]; n<- dim(Xfull)[1]
	mBF<- matrix(0, nc=k+3, nr=2^k, byr=T); colnames(mBF)<- c(colnames(Xfull),"dev","tBF","postprob")
	#first the null:
	l0<- logLik(survreg(Surv(y, delta)~1, dist="lognormal"))
	mBF[2^k, k+1]<- 0; mBF[2^k, k+2]<- 1
	for (i in 1:(2^k-1)){
		#binary expression of the model
		mBF[i, 1:k]<- integer.base.b(i, k)
		#obtain the deviance
		fit<- survreg(Surv(y, delta)~., data=as.data.frame(cbind(y, delta, Xfull[, mBF[i, 1:k]==1])), dist="lognormal")
		#deviance:
		mBF[i, k+1]<- 2*(logLik(fit)-l0)		
		dj<- sum(mBF[i, 1:k])
	
		#ZS adapted test-based BF
		a<- 1/2; b<- (sum(delta)+3)/2
		ap<- a+dj/2; bp<- b+mBF[i, k+1]/2
		
		lfact<- a*log(b)-lgamma(a)+pgamma(q=b, shape=a, scale=1)-
		        (ap*log(bp)-lgamma(ap)+pgamma(q=bp, shape=ap, scale=1))
				
		mBF[i, k+2]<- exp(lfact+mBF[i, k+1]/2)
	}
	  #posterior probs
		mBF[,k+3]<- mBF[, k+2]/((k+1)*choose(k, apply(mBF[,1:k], 1, sum)))
		mBF[,k+3]<- mBF[, k+3]/sum(mBF[,k+3])
	
	return(colSums(mBF[,1:k]*mBF[,k+3]))
}

