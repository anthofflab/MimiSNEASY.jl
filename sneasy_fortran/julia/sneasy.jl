const doeclimlib = joinpath(dirname(@__FILE__), "..", "sneasy", "doeclim")
const ccmlib = joinpath(dirname(@__FILE__), "..", "sneasy", "ccm")
const sneasylib = joinpath(dirname(@__FILE__), "..", "sneasy", "sneasy")

function run_fortran_doeclim(S::Float64, kappa::Float64, forcing::Vector{Float64})
	n = length(forcing)

	mod_time = zeros(n)
	mod_temp = zeros(n)
	mod_heatflux_mixed = zeros(n)
	mod_heatflux_interior = zeros(n)

	fout = ccall( (:run_doeclim_, doeclimlib),
		Int32, (Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		&n,
		mod_time,
		forcing,
		&S,
		&kappa,
		mod_temp,
		mod_heatflux_mixed,
		mod_heatflux_interior
		)

	return mod_temp, mod_heatflux_mixed, mod_heatflux_interior
end

function init_fortran_ccm()
	# load ocean anomaly table.
	oceanomtable = readdlm(joinpath(dirname(@__FILE__), "..", "sneasy", "anomtable.txt"))'
	ntab1, ntab2 = size(oceanomtable)

	ccall((:init_ccm_, ccmlib),Int,(Ptr{Int}, Ptr{Int}, Ptr{Float64}), &ntab1, &ntab2, oceanomtable)

	return
end

function fin_fortran_ccm()
	ccall((:fin_ccm_, ccmlib),Int,())
	return
end

function run_fortran_ccm(
		Climate_sens::Float64,
		Soil_respiration::Float64,
		Carbon_fertilization::Float64,
		Thermocline_diff::Float64,
		initial_co2::Float64,
		temp_forcing::Vector{Float64},
		CO2_emis_forcing::Vector{Float64})

	n = length(temp_forcing)

	atmco2_out = zeros(n)

	fout = ccall( (:run_ccm_, ccmlib),
		Int32, (Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		&n,
		temp_forcing,
		CO2_emis_forcing,
		&Climate_sens,
		&Soil_respiration,
		&Carbon_fertilization,
		&Thermocline_diff,
		&initial_co2,
		atmco2_out
		)

	return atmco2_out
end

function init_fortran_sneasy()
	# load ocean anomaly table.
	oceanomtable = readdlm(joinpath(dirname(@__FILE__), "..", "sneasy", "anomtable.txt"))'
	ntab1, ntab2 = size(oceanomtable)

	ccall((:init_sneasy_, sneasylib),Int,(Ptr{Int}, Ptr{Int}, Ptr{Float64}), &ntab1, &ntab2, oceanomtable)

	return
end

function fin_fortran_sneasy()
	ccall((:fin_sneasy_, sneasylib),Int,())
	return
end

function run_fortran_sneasy!(
	MOC_strength::Vector{Float64},
	radiative_forc::Vector{Float64},
	ATM_CO2::Vector{Float64},
	atm_oc_flux::Vector{Float64},
	GL_surface_temp::Vector{Float64},
	GL_ocean_heat::Vector{Float64},
	co2_emissions::Vector{Float64},
	rf_aerosol::Vector{Float64},
	rf_nonco2::Vector{Float64},
	S::Float64,
	kappa::Float64,
	alpha::Float64,
	Q10::Float64,
	beta::Float64,
	eta::Float64,
	hydsens::Float64,
	init_CO2::Float64=280.,
	init_MOC::Float64=20.)

	n = length(co2_emissions)

	if length(rf_aerosol)!=n || length(rf_nonco2)!=n
		error("All input vectors must have the same length")
	end

	tstep = 1.0     # years

	# call Fortran SNEASY
	ccall((:run_sneasy_model_, sneasylib), Int32,
		(Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		&n, &tstep, &S, &kappa, &alpha, &Q10, &beta, &eta, &hydsens, &init_MOC, &init_CO2, co2_emissions, rf_aerosol, rf_nonco2, MOC_strength, radiative_forc, ATM_CO2, atm_oc_flux, GL_surface_temp, GL_ocean_heat)

	return
end

function run_fortran_sneasy(
		co2_emissions::Vector{Float64},
		rf_aerosol::Vector{Float64},
		rf_nonco2::Vector{Float64},
		S::Float64, kappa::Float64,
		alpha::Float64,
		Q10::Float64,
		beta::Float64,
		eta::Float64,
		hydsens::Float64,
		init_CO2::Float64=280.,
		init_MOC::Float64=20.)

	n = length(co2_emissions)

	if length(rf_aerosol)!=n || length(rf_nonco2)!=n
		error("All input vectors must have the same length")
	end

	tstep = 1.0     # years

	# Allocate vectors for results
	MOC_strength = zeros(n)
	radiative_forc = zeros(n)
	ATM_CO2 = zeros(n)
	atm_oc_flux = zeros(n)
	GL_surface_temp = zeros(n)
	GL_ocean_heat = zeros(n)

	run_fortran_sneasy!(MOC_strength, radiative_forc, ATM_CO2, atm_oc_flux, GL_surface_temp, GL_ocean_heat, co2_emissions, rf_aerosol, rf_nonco2, S, kappa,	alpha,Q10, beta, eta, hydsens, init_CO2, init_MOC)

	return MOC_strength, radiative_forc, ATM_CO2, atm_oc_flux, GL_surface_temp, GL_ocean_heat
end
