using Base.Test
using Mimi
using DataFrames

include("../../sneasy/julia/sneasy.jl")
include("../src/doeclim.jl")

df = readtable("../../sneasy/data/forcing_rcp85.txt", separator=' ');
df = DataFrame(year=df[:year], rf=df[:co2]+df[:aerosol_direct]+df[:aerosol_indirect]+df[:ghg_nonco2]+df[:solar]+df[:volcanic]+df[:other]);
total_forcing = convert(Array, df[:rf]);

deltat = 1.0


function mimidoeclim(t2co, kappa, forcing)

    m = Model(;)

    setindex(m, :time, length(total_forcing))

    addcomponent(m, doeclimcomponent.doeclim, :doeclim)

	setparameter(m, :doeclim, :t2co, t2co)
    setparameter(m, :doeclim, :kappa, kappa)
    setparameter(m, :doeclim, :deltat, deltat)
    setparameter(m, :doeclim, :forcing, forcing)

    run(m)
    return m[:doeclim, :temp], m[:doeclim, :heatflux_mixed], m[:doeclim, :heatflux_interior]
end



#Note on Parameter/Variable Values
#	t2co				=	climate sensitivity to 2xCO2 (K)
#	kappa				=	vertical ocean diffusivity (cm^2 s^-1)
#	deltat				=	timestep
#	forcing				=	total radiative forcing
#	temp				=	global mean temperature anomaly (K), relative to preindustrial
#	heatflux_mixed		=	heat uptake of the mixed layer (W/m^2)
# 	heatflux_interior	=	heat uptake of the interior ocean (W/m^2)
#
#	format for run_fortran_doeclim = run_fortran_doeclim(t2co, Kappa, radiative_forcing)

#Run Mimi Verseion of doeclim
m_temp, m_heatflux_mixed, m_heatflux_interior = mimidoeclim(2.0, 1.1, total_forcing);

#Run fortran version of doeclim
f_temp, f_heatflux_mixed, f_heatflux_interior = run_fortran_doeclim(2.0, 1.1, total_forcing);


#Test Precision
Precision = 1.0e-6

@test_approx_eq_eps maxabs(m_temp .- f_temp) 0. Precision
@test_approx_eq_eps maxabs(m_heatflux_mixed .- f_heatflux_mixed) 0. Precision
@test_approx_eq_eps maxabs(m_heatflux_interior .- f_heatflux_interior) 0. Precision
