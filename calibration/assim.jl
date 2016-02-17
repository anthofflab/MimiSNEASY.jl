using Distributions
using DataFrames
using DataFramesMeta

include("loadobs.jl")
include(joinpath(dirname(@__FILE__), "..", "..", "sneasy", "julia", "sneasy.jl"))
include("run_mimi_sneasy.jl")

function sneasy_load_data()
	# load emissions data time series.
	df_emissions = readtable(joinpath(dirname(@__FILE__),"data/RCP85_EMISSIONS.csv"))
	rename!(df_emissions, :YEARS, :year)

	df_forcings = readtable(joinpath(dirname(@__FILE__),"data/forcing_rcp85.txt"), separator=' ')

	df = join(df_emissions, df_forcings, on=:year)

	df = DataFrame(year=df[:year], co2=df[:FossilCO2]+df[:OtherCO2], rf_aerosols=df[:aerosol_direct]+df[:aerosol_indirect], rf_other=df[:ghg_nonco2]+df[:solar]+df[:volcanic]+df[:other])

	return df
end

# log likelihood for a zero-mean AR1 process (innovation variance sigma^2, lag-1 autocorrelation coefficient rho1); handles initial value assuming process stationarity
function loglar1(r, σ, ρ1)
	n = length(r)
	σ_proc = σ/sqrt(1-ρ1^2) # stationary process variance sigma.proc^2

  	logl = logpdf(Normal(0, σ_proc), r[1])

  	for i=2:length(r)
		w = r[i] - ρ1 * r[i-1] # whitened residuals
		logl = logl + logpdf(Normal(0, σ), w)
	end
  	return logl
end

# (log) prior for model/statistical parameters
function log_pri(p)
	S = p[1]
	κ = p[2]
	α = p[3]
	Q10 = p[4]
	beta = p[5]
	eta = p[6]
	hydsens = p[7]
	T0 = p[8]
	H0 = p[9]
	CO20 = p[10]
	MOC0 = p[11]
	σ_temp = p[12]
	σ_ocheat = p[13]
	σ_co2inst = p[14]
	σ_co2ice = p[15]
	ρ_temp = p[16]
	ρ_ocheat = p[17]
	ρ_co2inst = p[18]

	prior_S1 = NormalInverseGaussian(1.7, 1.8, 1.2, 2.3)
	prior_S2 = NormalInverseGaussian(1.3, 1.9, 1.0, 3.3)
	prior_κ = LogNormal(1.1, 0.3)
	prior_α = TriangularDist(0., 3., 1.)
	prior_Q10 = Uniform(0., 5)
	prior_beta = Uniform(0., 1)
	prior_eta = Uniform(0., 200)
	prior_hydsens = Uniform(0., 0.06)
	prior_T0 = Normal()
	prior_H0 = Uniform(-100, 0)
	prior_CO20 = Uniform(280, 295)
	prior_MOC0 = Uniform(10, 30)
	prior_σ_temp = Uniform(0, 0.2)
	prior_σ_ocheat = Uniform(0, 4)
	prior_σ_co2inst = Uniform(0, 1)
	prior_σ_co2ice = Uniform(0, 10)
	prior_ρ_temp = Uniform(0, 0.99)
	prior_ρ_ocheat = Uniform(0, 0.99)
	prior_ρ_co2inst = Uniform(0, 0.99)

	# For climate sensitivyt the R code uses a truncated distribution in the R
	# code with support 0 to Inf. The best way to fix this would be to add
	# support for a truncated NormalInverseGaussian to the Distributions package.
	# For now, the code below just checks the bounds manually, and then uses
	# the pdf of the NormalInverseGaussian without any truncation. Technically
	# that is not the correct pdf, but for the MCMC algorithm that doesn't matter.
	lpri = -Inf
	if S>0.
        lpri = logpdf(prior_S1, S) + logpdf(prior_S2, S) + logpdf(prior_κ, κ) + logpdf(prior_α, α) + logpdf(prior_Q10, Q10) + logpdf(prior_beta, beta) + logpdf(prior_eta, eta) + logpdf(prior_hydsens, hydsens) + logpdf(prior_T0, T0) + logpdf(prior_H0, H0) + logpdf(prior_CO20, CO20) + logpdf(prior_MOC0, MOC0) + logpdf(prior_σ_temp, σ_temp) + logpdf(prior_σ_ocheat, σ_ocheat) + logpdf(prior_σ_co2inst, σ_co2inst) + logpdf(prior_σ_co2ice, σ_co2ice) + logpdf(prior_ρ_temp, ρ_temp) + logpdf(prior_ρ_ocheat, ρ_ocheat) + logpdf(prior_ρ_co2inst, ρ_co2inst)
	end

	# This is just to emulate the R code right now, i.e. if the prior is finite
	# we are recalculating the prior to be just some parameters, see issue #9
	# TODO remove this whole if clause once things are cross checked
	if isfinite(lpri)
		lpri = logpdf(prior_S1, S) + logpdf(prior_S2, S) + logpdf(prior_κ, κ) + logpdf(prior_α, α) + logpdf(prior_T0, T0)
	end

	return lpri
end

function construct_log_post(f_run_model, endyear=2010; assim_temp=true, assim_ocheat=true, assim_co2inst=true, assim_co2ice=true, assim_ocflux=true)
	df_forcing = sneasy_load_data()

	df_obs = loaddata()
	df = join(df_forcing, df_obs, on=:year, kind=:outer)
	sort!(df, cols=[:year])

	df = @where(df, :year .<= 2010)

	f_co2 = convert(Array, df[:co2])
	f_aerosols = convert(Array, df[:rf_aerosols])
	f_other = convert(Array, df[:rf_other])

	obs_temperature = df[:obs_temperature]
  	obs_temperature_indiceswithdata = Array(Int,0)
	for i=1:length(obs_temperature)
		if !isna(obs_temperature[i])
			push!(obs_temperature_indiceswithdata, i)
		end
	end
  	tempvar_temperature_res = zeros(length(obs_temperature_indiceswithdata))

	obs_ocheat = df[:obs_ocheat]
  	obs_ocheat_indiceswithdata = Array(Int,0)
	for i=1:length(obs_ocheat)
		if !isna(obs_ocheat[i])
			push!(obs_ocheat_indiceswithdata, i)
		end
	end
  	tempvar_ocheat_res = zeros(length(obs_ocheat_indiceswithdata))

	obs_co2inst = df[:obs_co2inst]
  	obs_co2inst_indiceswithdata = Array(Int,0)
	for i=1:length(obs_co2inst)
		if !isna(obs_co2inst[i])
			push!(obs_co2inst_indiceswithdata, i)
		end
	end
  	tempvar_co2inst_res = zeros(length(obs_co2inst_indiceswithdata))

    obs_co2ice = df[:obs_co2ice]
    obs_co2ice_indiceswithdata = Array(Int,0)
    for i=1:length(obs_co2ice)
        if !isna(obs_co2ice[i])
            push!(obs_co2ice_indiceswithdata, i)
        end
    end
    mean_co2ice = zeros(length(obs_co2ice_indiceswithdata))

    obs_ocflux = df[:obs_ocflux]
    obs_ocflux_err = df[:obs_ocflux_err]
    obs_ocflux_indiceswithdata = Array(Int,0)
    for i=1:length(obs_ocflux)
        if !isna(obs_ocflux[i])
            push!(obs_ocflux_indiceswithdata, i)
        end
    end

	init_fortran_sneasy()

	n = length(f_co2)

	# Allocate vectors for results
	MOC_strength = zeros(n)
	radiative_forc = zeros(n)
	model_co2 = zeros(n)
	atm_oc_flux = zeros(n)
	model_temperature = zeros(n)
	model_ocheat = zeros(n)

	# (log) likelihood for observations, assuming residual independence between data sets
	function sneays_log_lik(p)
		S = p[1]
		κ = p[2]
		α = p[3]
		Q10 = p[4]
		beta = p[5]
		eta = p[6]
		hydsens = p[7]
		T0 = p[8]
		H0 = p[9]
		CO20 = p[10]
		MOC0 = p[11]
		σ_temp = p[12]
		σ_ocheat = p[13]
		σ_co2inst = p[14]
		σ_co2ice = p[15]
		ρ_temp = p[16]
		ρ_ocheat = p[17]
		ρ_co2inst = p[18]

		f_run_model(MOC_strength, radiative_forc, model_co2, atm_oc_flux, model_temperature, model_ocheat,
			f_co2, f_aerosols, f_other,
			S, κ, α, Q10, beta, eta, hydsens, CO20, MOC0)

		llik_temp = 0.
		if assim_temp
			for (i, index)=enumerate(obs_temperature_indiceswithdata)
				tempvar_temperature_res[i] = obs_temperature[index] - (model_temperature[index] + T0)
			end
			llik_temp = loglar1(tempvar_temperature_res, σ_temp, ρ_temp)
		end

		llik_ocheat = 0.
		if assim_ocheat
			for (i, index)=enumerate(obs_ocheat_indiceswithdata)
				tempvar_ocheat_res[i] = obs_ocheat[index] - (model_ocheat[index] + H0)
			end
			llik_ocheat = loglar1(tempvar_ocheat_res, σ_ocheat, ρ_ocheat)
		end

		llik_co2inst = 0.
		if assim_co2inst
			for (i, index)=enumerate(obs_co2inst_indiceswithdata)
				tempvar_co2inst_res[i] = obs_co2inst[index] - model_co2[index]
			end
			llik_co2inst = loglar1(tempvar_co2inst_res, σ_co2inst, ρ_co2inst)
		end

        llik_co2ice = 0.
        if assim_co2ice
            for (i, index)=enumerate(obs_co2ice_indiceswithdata)
                mean_co2ice[i] = mean(model_co2[index + (-4:3)])
                llik_co2ice = llik_co2ice + logpdf(Normal(mean_co2ice[i], σ_co2ice), obs_co2ice[index])
            end
        end

        llik_ocflux = 0.
        if assim_ocflux
            for (i,index) = enumerate(obs_ocflux_indiceswithdata)
                llik_ocflux =  llik_ocflux + logpdf(Normal(atm_oc_flux[index], obs_ocflux_err[index]), obs_ocflux[index])
            end
        end

		#llik.moc = 0
		#if(assim.moc & !is.null(oidx.moc))
		#	llik.moc = sum(dnorm(obs.moc[oidx.moc], model.out$moc[midx.moc], obs.moc.err[oidx.moc], log=TRUE))

		#llik = llik.temp + llik.ocheat + llik.co2inst + llik.co2inst + llik.co2ice + llik.ocflux + llik.moc # assume the residual time series are independent

		llik = llik_temp + llik_ocheat + llik_co2inst + llik_co2ice + llik_ocflux

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
function proposal_matrix(chain; mult=1., beta=0.05)
#	# mult = overall scale factor to adjust all step sizes
#	# beta = relative influence of nugget term
	prechain = chain.value'
#
	p = size(prechain, 2)
	precov = cov(prechain)

	propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*eye(p)/p
	propcov = mult*propcov

	mat = chol(propcov)'

	return convert(Array, mat)
end
