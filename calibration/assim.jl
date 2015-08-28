using Distributions
using DataFrames
using DataFramesMeta
using Mimi

include("loadobs.jl")
include("../src/sneasy.jl")

# log likelihood for a zero-mean AR1 process (innovation variance sigma^2, lag-1 autocorrelation coefficient rho1); handles initial value assuming process stationarity
function loglar1(r, σ, ρ1)
	n = length(r)
	σ_proc = σ/sqrt(1-ρ1^2) # stationary process variance sigma.proc^2

  nd = Normal(0., σ_proc)

  logl = logpdf(nd, r[1])

  for i=2:length(r)
		w = r[i] - ρ1 * r[i-1] # whitened residuals
		logl = logl + logpdf(nd, w)
	end

  return logl
end

# (log) prior for model/statistical parameters
function log_pri(p)
	S = p[1]
	#κ = p[2]
	#α = p[3]
	#Q10 = p[4]
	#beta = p[5]
	#eta = p[6]
	#hydsens = p[7]
	#T0 = p[8]
	#H0 = p[9]
	#CO20 = p[10]
	#MOC0 = p[11]
	σ_temp = p[2]
	#σ_ocheat = p[13]
	#σ_co2inst = p[14]
	#σ_co2ice = p[15]
	ρ_temp = p[3]
	#ρ_ocheat = p[17]
	#ρ_co2inst = p[18]

	#bound_lower = [0,0,0,0,0,0,0,-Inf,-100,280,10,0,0,0,0,0,0,0]
	#bound_upper = [Inf,Inf,3,5,1,200,0.06,Inf,0,295,30,0.2,4,1,10,0.99,0.99,0.99]
	#in.range = all(p > bound.lower) & all(p < bound.upper)
  prior_S = Normal(3., 2.)

	lpri_cs = insupport(prior_S, S) ? logpdf(prior_S, S) : -Inf

	lpri = lpri_cs

#	if(in.range) {
#		lpri.cs = dnig(S, alpha=1.8, delta=2.3, beta=1.2, mu=1.7, log=TRUE) + dnig(S, alpha=1.9, delta=3.3, beta=1.0, mu=1.3, log=TRUE) # NIG mean state + LGM prior
#		lpri.kappa = dlnorm(kappa, meanlog=1.1, sdlog=0.3, log=TRUE) # lognormal prior determined from UVic biogeochemical tracers
#		if(alpha < 1) lpri.alpha = log(2/3*alpha) else lpri.alpha = log(1-1/3*alpha) # triangular prior approximately based on IPCC AR4 WG1 Fig SPM.2
#		lpri.T0 = dnorm(T0, 0, 1, log=TRUE) # normal prior on T0
#		lpri = lpri.cs + lpri.kappa + lpri.alpha + lpri.T0
#	} else {
#		lpri = -Inf
#	}
#
#	names(lpri) = NULL
	lpri
end

function construct_loglikelihood(endyear=2010)
	df = loaddata()

	df = @where(df, :year .<= 2010)

  forcing_emis_co2 = convert(Array, df[:forcing_emis_co2])
  forcing_rf_nonco2 = convert(Array, df[:forcing_rf_nonco2])

	obs_temperature = df[:obs_temperature]
  obs_temperature_indiceswithdata = Array(Int,0)
	for i=1:length(obs_temperature)
		if !isna(obs_temperature[i])
			push!(obs_temperature_indiceswithdata, i)
		end
	end
  tempvar_temperature_res = zeros(length(obs_temperature_indiceswithdata))

	sneasy_model = getsneasy()

	model_temperature = sneasy_model[:doeclim, :temp]

	# (log) likelihood for observations, assuming residual independence between data sets
	function sneays_log_lik(p)
		S = p[1]
  	#κ = p[2]
  	#α = p[3]
		#Q10 = p[4]
  	#beta = p[5]
  	#eta = p[6]
		#hydsens = p[7]
		#T0 = p[8]
  	#H0 = p[9]
  	#CO20 = p[10]
  	#MOC0 = p[11]
		σ_temp = p[2]
  	#σ_ocheat = p[13]
  	#σ_co2inst = p[14]
  	#σ_co2ice = p[15]
		ρ_temp = p[3]
  	#ρ_ocheat = p[17]
  	#ρ_co2inst = p[18]

		setparameter(sneasy_model, :ccm, :Clim_sens, S)

		run(sneasy_model)

		llik_temp = 0.
		#if(assim.temp & !is.null(oidx.temp)) {
		for (i, index)=enumerate(obs_temperature_indiceswithdata)
			tempvar_temperature_res[i] = obs_temperature[index] - model_temperature[index]
		end
		lik_temp = loglar1(tempvar_temperature_res, σ_temp, ρ_temp)
#			resid_temp = obs.temp[oidx.temp] - (model.out$temp[midx.temp]+T0)
#			llik.temp = logl.ar1(resid.temp, sigma.temp, rho.temp) # AR(1)
#		}

	#=	llik.ocheat = 0
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

	llik = llik.temp + llik.ocheat + llik.co2inst + llik.co2inst + llik.co2ice + llik.ocflux + llik.moc # assume the residual time series are independent
=#
		llik = llik_temp

		return llik
	end

	# (log) posterior distribution:  posterior ~ likelihood * prior
	function log_post(p)
		lpri = log_pri(p)

		lpost = isfinite(lpri) ? sneays_log_lik(p) + lpri : -Inf

		return lpost
	end

	return log_post
end




# calculate a scale matrix to transform a iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
#proposal.matrix = function(prechain, mult=1, beta=0.05)
#{
#	# mult = overall scale factor to adjust all step sizes
#	# beta = relative influence of nugget term
#
#	p = ncol(prechain)
#	precov = cov(prechain)

#	propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p
#	propcov = mult*propcov

#	mat = t(chol(propcov))
#}
