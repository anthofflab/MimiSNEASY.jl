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


# load forcing time series
forcing = read.table("../sneasy/forcing.txt", col.names=c("year","co2","nonco2.land","nonco2.ocean","aerosol.land","aerosol.ocean","solar.land","solar.ocean","volc.land","volc.ocean","tot.land","tot.ocean"))
mod.time = forcing$year
endyear = mod.time[length(mod.time)]

# calculate total radiative forcing from individual land/ocean forcings
total.forcing = function(forcing, alpha)
{
	flnd = 0.29 # area land fraction
	
	forcing.land = forcing$co2 + forcing$nonco2.land + alpha*forcing$aerosol.land + forcing$solar.land + forcing$volc.land
	forcing.ocean = forcing$co2 + forcing$nonco2.ocean + alpha*forcing$aerosol.ocean + forcing$solar.ocean + forcing$volc.ocean
	
	forcing.total = flnd*forcing.land + (1-flnd)*forcing.ocean
	
	return(forcing.total)
}

# convert annual ocean heat flux (W/m^2) to cumulative ocean heat content anomaly (10^22 J)
flux.to.heat = function(heatflux.mixed, heatflux.interior)
{
	flnd = 0.29 # area land fraction
	fso = 0.95 # ocean area fraction of interior 
	secs.per.year = 31556926
	earth.area = 510065600 * 10^6
	ocean.area = (1-flnd)*earth.area
	powtoheat = ocean.area*secs.per.year / 10^22 # in 10^22 J/yr
	
	heat.mixed = cumsum(heatflux.mixed) * powtoheat
	heat.interior = fso * cumsum(heatflux.interior) * powtoheat
	ocean.heat = heat.mixed + heat.interior
	
	return(list(ocean.heat=ocean.heat, heat.mixed=heat.mixed, heat.interior=heat.interior))
}




#####      load DOECLIM model shared library     #####
dyn.load("../sneasy/doeclim.dll")


# DOECLIM climate model (0D EBM atmosphere + 1D diffusive ocean)
# inputs: climate sensitivity (S), ocean vertical diffusivity (kappa), aerosol forcing scale factor (alpha)
# outputs: annual global temperature (temp, K) and total ocean heat (ocheat, 10^22 J), mixed-layer and interior ocean heat (ocheat.mixed and ocheat.interior, 10^22 J), atmosphere->mixed and mixed->interior heat fluxes (ocheatflux.mixed and ocheatflux.interior, W/m^2)
doeclim = function(S, kappa, alpha)
{	
	n = length(mod.time)
	forcing.total = total.forcing(forcing, alpha)
		
	# call Fortran DOECLIM
	# doeclim.so must be already dynamically loaded (see above this function)
	fout = .Fortran( "run_doeclim",
			ns = n,
			time_out = as.double(mod.time),
			forcing_in = as.double(forcing.total),
			t2co_in = as.double(S),
			kappa_in = as.double(kappa),
			temp_out = as.double(rep(0,n)),
			heatflux_mixed_out = as.double(rep(0,n)),
			heatflux_interior_out = as.double(rep(0,n))
		)
	
	ocheat = flux.to.heat(fout$heatflux_mixed, fout$heatflux_interior)

	model.output = list(time=mod.time, temp=fout$temp, ocheat=ocheat$ocean.heat, ocheat.mixed=ocheat$heat.mixed, ocheat.interior=ocheat$heat.interior, ocheatflux.mixed = fout$heatflux_mixed, ocheatflux.interior = fout$heatflux_interior)
	
	return(model.output)
}
