using Base.Test

include("../src/sneasy.jl")
include("../../sneasy/julia/sneasy.jl")

#DATA
#Note: f_nonco2forcing can be passed into mimisneasy, but the fortran version requires the forcing inputs to be split up.
df = readtable("../../sneasy/data/RCP85_EMISSIONS.csv");
f_co2emissions = convert(Array, (df[:FossilCO2]+df[:OtherCO2]));
f_nonco2forcing = readtable("../../sneasy/data/forcing_rcp85.txt", separator = ' ', header=true);
f_aerosols=convert(Array, (f_nonco2forcing[:aerosol_direct]+f_nonco2forcing[:aerosol_indirect]));
f_other = convert(Array, (f_nonco2forcing[:ghg_nonco2]+f_nonco2forcing[:solar]+f_nonco2forcing[:volcanic]+f_nonco2forcing[:other]));



function mimisneasy(f_co2emissions, f_nonco2forcing, t2co, kappa, alpha, Q10, Beta, Eta)
    m = Model(;)

    setindex(m, :time, nrow(f_nonco2forcing))

    addcomponent(m, radforccomponent.radforc, :radforc)
    addcomponent(m, doeclimcomponent.doeclim, :doeclim)
    addcomponent(m, ccmcomponent.ccm, :ccm)

    f_anomtable = readdlm("../../sneasy/sneasy/anomtable.txt");

    # Timesteps
    deltat = 1.0
    MOC = "Not calculated for Mimi"

    setparameter(m, :doeclim, :t2co, t2co)
    setparameter(m, :doeclim, :kappa, kappa)
    setparameter(m, :doeclim, :deltat, deltat)

    setparameter(m, :ccm, :deltat, deltat)
    setparameter(m, :ccm, :Q10, Q10)
    setparameter(m, :ccm, :Beta, Beta)
    setparameter(m, :ccm, :Eta, Eta)
    setparameter(m, :ccm, :atmco20, 280.) #Set to 280 in run_fortran_sneasy
    setparameter(m, :ccm, :CO2_emissions, f_co2emissions)
    setparameter(m, :ccm, :anomtable, transpose(f_anomtable))

    setparameter(m, :radforc, :forcing_other, Vector(f_nonco2forcing[:other]))
    setparameter(m, :radforc, :forcing_volcanic, Vector(f_nonco2forcing[:volcanic]))
    setparameter(m, :radforc, :forcing_solar, Vector(f_nonco2forcing[:solar]))
    setparameter(m, :radforc, :forcing_ghg_nonco2, Vector(f_nonco2forcing[:ghg_nonco2]))
    setparameter(m, :radforc, :forcing_aerosol_direct, Vector(f_nonco2forcing[:aerosol_direct]))
    setparameter(m, :radforc, :forcing_aerosol_indirect, Vector(f_nonco2forcing[:aerosol_indirect]))
    setparameter(m, :radforc, :alpha, alpha)
    setparameter(m, :radforc, :deltat, deltat)

    connectparameter(m, :doeclim, :forcing, :radforc, :rf)
    connectparameter(m, :ccm, :temp, :doeclim, :temp)
    connectparameter(m, :radforc, :atmco2, :ccm, :atmco2)

    run(m)
    return MOC, m[:radforc, :rf], m[:ccm, :atmco2], m[:ccm, :atm_oc_flux], m[:doeclim, :temp], m[:doeclim, :heat_interior]
end

#Note on Parameter Values
#	t2co	=	climate sensitivity to 2xCO2 (K)
#	kappa	=	vertical ocean diffusivity (cm^2 s^-1)
#	Q10		=	Respiration Temperature sens.
#	Beta	=	Carbon Fertilization param.
#	Eta		=	Thermocline transfer velocity
#	alpha	=	aerosol amplification
#
#format for run_fortran_sneasy = run_fortran_sneasy(co2_emissions_forcing, aerosol_forcing, other_forcing, t2co, kappa, alpha, Q10, Beta, Eta)

m_MOC, m_radforc, m_atmco2, m_atmocflux, m_surfacetemp, m_heatinterior = mimisneasy(f_co2emissions, f_nonco2forcing, 2.0, 1.1, 1., 1.311, 0.502, 17.7);

#Note, last parameter (0.047) in run_fortran_sneasy is for the MOC and does not influence results being compared.
init_fortran_sneasy()
f_MOC, f_radforc, f_atmco2, f_atmocflux, f_surfacetemp, f_heatinterior = run_fortran_sneasy(f_co2emissions, f_aerosols, f_other, 2.0, 1.1, 1., 1.311, 0.502, 17.7, 0.047);
fin_fortran_sneasy()

#Test Precision
Precision = 1.0e-3

@test_approx_eq_eps maxabs(m_radforc .- f_radforc) 0. Precision
@test_approx_eq_eps maxabs(m_atmco2 .- f_atmco2) 0. Precision
@test_approx_eq_eps maxabs(m_atmocflux .- f_atmocflux) 0. Precision
@test_approx_eq_eps maxabs(m_surfacetemp .- f_surfacetemp) 0. Precision
@test_approx_eq_eps maxabs(m_heatinterior .- f_heatinterior) 0. Precision

