# doeclim.R
#
# Nathan M. Urban (nurban@psu.edu)
# Department of Geosciences, Penn State
#
# Implements DOECLIM, a simple climate model
#
# DOECLIM is a 0-dimensional energy balance model (EBM) for the
# atmosphere coupled to a 1-dimensional diffusive ocean.  The
# model outputs temperature and ocean heat content time series
# as well as ocean heat fluxes.  See:
#
#  Elmar Kriegler, "Imprecise probability analysis for integrated
#  assessment of climate change", Ph.D. thesis, Potsdam (2005).
#   http://opus.kobv.de/ubp/volltexte/2005/561/
#  (Fortran port by Marlos Goes and Nathan Urban.)
#
# The model is implemented in Fortran and called from R.  This
# file also contains functions to load and process forcing data
# and model output.  Any pre/post-processing of input/output
# should be done in R or otherwise separate from the main
# model subroutine, which for computational efficiency
# should not perform any file I/O itself.  The Fortran model
# must be implemented as a standalone subroutine.
#
# For further information on R/Fortran coupling, see:
#
#   http://www.stat.umn.edu/~charlie/rc/
#   http://math.acadiau.ca/ACMMaC/howtos/Fortran_R.html
#   http://cran.r-project.org/doc/manuals/R-exts.pdf (chapter 5)
using DataFrames
using Base.Test
include("../src/sneasy.jl")

# load forcing time series
forcing = readtable("forcing.txt", separator=' ', names=[:year,:co2,:nonco2_land,:nonco2_ocean,:aerosol_land,:aerosol_ocean,:solar_land,:solar_ocean,:volc_land,:volc_ocean,:tot_land,:tot_ocean], header=false)
forcing_time = forcing[:year]
endyear = forcing_time[end]

# calculate total radiative forcing from individual land/ocean forcings
function total_forcing(forcing, alpha)
	flnd = 0.29 # area land fraction

	forcing_land = forcing[:co2] + forcing[:nonco2_land] + alpha*forcing[:aerosol_land] + forcing[:solar_land] + forcing[:volc_land]
	forcing_ocean = forcing[:co2] + forcing[:nonco2_ocean] + alpha*forcing[:aerosol_ocean] + forcing[:solar_ocean] + forcing[:volc_ocean]

	forcing_total = flnd .* forcing_land + (1-flnd) .* forcing_ocean

	return forcing_total
end

# convert annual ocean heat flux (W/m^2) to cumulative ocean heat content anomaly (10^22 J)
#flux.to.heat = function(heatflux.mixed, heatflux.interior)
#{
#	flnd = 0.29 # area land fraction
#	fso = 0.95 # ocean area fraction of interior
#	secs.per.year = 31556926
#	earth.area = 510065600 * 10^6
#	ocean.area = (1-flnd)*earth.area
#	powtoheat = ocean.area*secs.per.year / 10^22 # in 10^22 J/yr
#
#	heat.mixed = cumsum(heatflux.mixed) * powtoheat
#	heat.interior = fso * cumsum(heatflux.interior) * powtoheat
#	ocean.heat = heat.mixed + heat.interior
#
#	return(list(ocean.heat=ocean.heat, heat.mixed=heat.mixed, heat.interior=heat.interior))
#}

# DOECLIM climate model (0D EBM atmosphere + 1D diffusive ocean)
# inputs: climate sensitivity (S), ocean vertical diffusivity (kappa), aerosol forcing scale factor (alpha)
# outputs: annual global temperature (temp, K) and total ocean heat (ocheat, 10^22 J), mixed-layer and interior ocean heat (ocheat.mixed and ocheat.interior, 10^22 J), atmosphere->mixed and mixed->interior heat fluxes (ocheatflux.mixed and ocheatflux.interior, W/m^2)
function doeclim(S, kappa, alpha)
	n = length(forcing_time)
	forcing_total = total_forcing(forcing, alpha)

	mod_time = zeros(n)
	mod_temp = zeros(n)
	mod_heatflux_mixed = zeros(n)
	mod_heatflux_interior = zeros(n)

	# call Fortran DOECLIM
	fout = ccall( (:run_doeclim_, "doeclim"),
		Int32, (Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		&n,
		mod_time,
		convert(Array, forcing_total),
		&S,
		&kappa,
		mod_temp,
		mod_heatflux_mixed,
		mod_heatflux_interior
		)

	#ocheat = flux.to.heat(fout$heatflux_mixed, fout$heatflux_interior)

	#model.output = list(time=mod.time, temp=fout$temp, ocheat=ocheat$ocean.heat, ocheat.mixed=ocheat$heat.mixed, ocheat.interior=ocheat$heat.interior, ocheatflux.mixed = fout$heatflux_mixed, ocheatflux.interior = fout$heatflux_interior)

	#return(model.output)
	return mod_temp, mod_heatflux_mixed, mod_heatflux_interior
end

function juliaversion(S, kappa, alpha)
	forcing_total = total_forcing(forcing, alpha)

	m = Model()

	setindex(m, :time, length(forcing_time))

	addcomponent(m, doeclimcomponent.doeclim)
	setparameter(m, :doeclim, :t2co, S)
	setparameter(m, :doeclim, :kappa, kappa)
	setparameter(m, :doeclim, :deltat, 1.)
	setparameter(m, :doeclim, :forcing, convert(Array, forcing_total))

	run(m)

	return m[:doeclim, :temp], m[:doeclim, :heatflux_mixed], m[:doeclim, :heatflux_interior]
end

temp_julia, mixed_julia, interior_julia = juliaversion(2., 1.1, 0.6)
temp_f, mixed_f, interior_f = doeclim(2., 1.1, 0.6)

@test maxabs(temp_julia-temp_f) < 0.0001
