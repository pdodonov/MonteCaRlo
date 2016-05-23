bootstrapT <- function (x, n=10000, alpha=0.05) { # requires a numeric vector as input
	#Critical values are calculated using bootstrap-t (Manly 1998, p. 56)
	Teta.boot=numeric(0)
	SEM.boot=numeric(0)
	Tb=numeric(0) 
	Teta.est=mean(x)
	SEM.est=sd(x)/sqrt(length(x))
	for (i in 1:n) {
		x.boot=sample(x,size=length(x),replace=TRUE,)
		Tb[i]=(mean(x.boot)-Teta.est)/(sd(x.boot)/sqrt(length(x)))
		}
	Tb=Tb[abs(Tb) != Inf]
	hist(Tb)
	#Calculating the alpha critical T values
	Tinf=quantile(Tb,alpha/2,na.rm=TRUE)
	Tsup=quantile(Tb,1-alpha/2,na.rm=TRUE)
	lim.inf=Teta.est-Tsup*SEM.est
	lim.sup=Teta.est-Tinf*SEM.est
	c(lim.inf,lim.sup)
	}

