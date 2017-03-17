# normal inverse Gaussian density
dnig = function(x, alpha=1.5, delta=1, beta=0, mu=0, log=FALSE)
{
	ret = log(alpha*delta) + log(besselK(alpha*sqrt(delta^2+(x-mu)^2),1)) - log(pi*sqrt(delta^2+(x-mu)^2)) + delta*sqrt(alpha^2-beta^2) + beta*(x-mu)
	if(!log) ret = exp(ret)
	return(ret)
}

# log likelihood for a zero-mean AR1 process (innovation variance sigma^2, lag-1 autocorrelation coefficient rho1); handles initial value assuming process stationarity
logl.ar1 = function(r,sigma,rho1)
{
	n = length(r)
	sigma.proc = sigma/sqrt(1-rho1^2) # stationary process variance sigma.proc^2

	logl = dnorm(r[1],sd=sigma.proc,log=TRUE)
	if(n>1) {
		w = r[2:n] - rho1*r[1:(n-1)] # whitened residuals
		logl = logl + sum(dnorm(w,sd=sigma,log=TRUE))
	}
}

# (log) likelihood for observations, assuming residual independence between data sets
log.lik = function(p) 
{
	S = p[1]; kappa = p[2]; alpha = p[3]
	Q10 = p[4]; beta = p[5]; eta = p[6]
	hydsens = p[7]
	T0 = p[8]; H0 = p[9]; CO20 = p[10]; MOC0 = p[11]
	sigma.temp = p[12]; sigma.ocheat = p[13]; sigma.co2inst = p[14]; sigma.co2ice = p[15]
	rho.temp = p[16]; rho.ocheat = p[17]; rho.co2inst = p[18]

	model.out = sneasy(S, kappa, alpha, Q10, beta, eta, hydsens, CO20, MOC0, endyear)
	
	llik.temp = 0
	if(assim.temp & !is.null(oidx.temp)) {
		resid.temp = obs.temp[oidx.temp] - (model.out$temp[midx.temp]+T0)
		llik.temp = logl.ar1(resid.temp, sigma.temp, rho.temp) # AR(1)
	}
	
	llik.ocheat = 0
	if(assim.ocheat & !is.null(oidx.ocheat)) {
		resid.ocheat = obs.ocheat[oidx.ocheat] - (model.out$ocheat[midx.ocheat]+H0)
		llik.ocheat = logl.ar1(resid.ocheat, sigma.ocheat, rho.ocheat) # AR(1)
	}
	
	llik.co2inst = 0
	if(assim.co2inst & !is.null(oidx.co2inst)) {
		resid.co2inst = obs.co2inst[oidx.co2inst] - model.out$co2[midx.co2inst]
		llik.co2inst = logl.ar1(resid.co2inst, sigma.co2inst, rho.co2inst) # AR(1)
	}
	
	llik.co2ice = 0
	if(assim.co2ice & !is.null(oidx.co2ice)) {
		mod.co2ice = rep(NA,length(midx.co2ice))
		for(i in 1:length(midx.co2ice))
			mod.co2ice[i] = mean(model.out$co2[midx.co2ice[i] + -4:3]) # mean of 8 years centered on ice core observation
		llik.co2ice = sum(dnorm(obs.co2ice[oidx.co2ice], mod.co2ice, sigma.co2ice, log=TRUE))
	}
	
	llik.ocflux = 0
	if(assim.ocflux & !is.null(oidx.ocflux))
		llik.ocflux = sum(dnorm(obs.ocflux[oidx.ocflux], model.out$ocflux[midx.ocflux], obs.ocflux.err[oidx.ocflux], log=TRUE))

	llik.moc = 0
	if(assim.moc & !is.null(oidx.moc))
		llik.moc = sum(dnorm(obs.moc[oidx.moc], model.out$moc[midx.moc], obs.moc.err[oidx.moc], log=TRUE))
	
	llik = llik.temp + llik.ocheat + llik.co2inst + llik.co2ice + llik.ocflux + llik.moc # assume the residual time series are independent
	llik
}

# (log) prior for model/statistical parameters
log.pri = function(p)
{
	S = p[1]; kappa = p[2]; alpha = p[3]
	Q10 = p[4]; beta = p[5]; eta = p[6]
	hydsens = p[7]
	T0 = p[8]; H0 = p[9]; CO20 = p[10]; MOC0 = p[11]
	sigma.temp = p[12]; sigma.ocheat = p[13]; sigma.co2inst = p[14]; sigma.co2ice = p[15]
	rho.temp = p[16]; rho.ocheat = p[17]; rho.co2inst = p[18]	
	bound.lower = c(0,0,0,0,0,0,0,-Inf,-100,280,10,0,0,0,0,0,0,0)
	bound.upper = c(Inf,Inf,3,5,1,200,0.06,Inf,0,295,30,0.2,4,1,10,0.99,0.99,0.99)
	in.range = all(p > bound.lower) & all(p < bound.upper)
		
	if(in.range) {
		lpri.cs = dnig(S, alpha=1.8, delta=2.3, beta=1.2, mu=1.7, log=TRUE) + dnig(S, alpha=1.9, delta=3.3, beta=1.0, mu=1.3, log=TRUE) # NIG mean state + LGM prior
		lpri.kappa = dlnorm(kappa, meanlog=1.1, sdlog=0.3, log=TRUE) # lognormal prior determined from UVic biogeochemical tracers
		if(alpha < 1) lpri.alpha = log(2/3*alpha) else lpri.alpha = log(1-1/3*alpha) # triangular prior approximately based on IPCC AR4 WG1 Fig SPM.2
		lpri.T0 = dnorm(T0, 0, 1, log=TRUE) # normal prior on T0
		lpri = lpri.cs + lpri.kappa + lpri.alpha + lpri.T0
	} else {
		lpri = -Inf
	}
	
	names(lpri) = NULL
	lpri
}

# (log) posterior distribution:  posterior ~ likelihood * prior
log.post = function(p)
{	
	lpri = log.pri(p)
	if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
		lpost = log.lik(p) + lpri
	} else {
		lpost = -Inf	
	}
	lpost
}

# calculate a scale matrix to transform a iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
proposal.matrix = function(prechain, mult=1, beta=0.05)
{
	# mult = overall scale factor to adjust all step sizes
	# beta = relative influence of nugget term

	p = ncol(prechain)
	precov = cov(prechain)
	
	propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p
	propcov = mult*propcov
	
	mat = t(chol(propcov))
}