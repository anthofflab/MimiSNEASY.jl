using Base.Test
using Mimi
using DataFrames
using Distributions

include("../../sneasy/julia/sneasy.jl")
include("../src/ccm.jl")

df = readtable("../../sneasy/data/RCP85_EMISSIONS.csv");
co2_forcing = convert(Array, (df[:FossilCO2] + df[:OtherCO2]));

#Load Ocean anomaly table (need to use transpose)
anomtable = readdlm("../../sneasy/sneasy/anomtable.txt");

#Create a temperature forcing file with the appropriate length
srand(123);
temp_forcing = Float64[0.8 * (1+0.0025)^t + rand(Normal(0., 0.2), 1)[1] for t = 1:length(co2_forcing)];

deltat = 1.0


function mimiccm(Q10, Beta, Eta, atmco20, temp_forcing, co2_forcing)

    m = Model()

    setindex(m, :time, length(co2_forcing))

 	addcomponent(m, ccmcomponent.ccm, :ccm)

    setparameter(m, :ccm, :deltat, deltat)
    setparameter(m, :ccm, :Q10, Q10)
    setparameter(m, :ccm, :Beta, Beta)
    setparameter(m, :ccm, :Eta, Eta)
    setparameter(m, :ccm, :atmco20, atmco20)
    setparameter(m, :ccm, :CO2_emissions, co2_forcing)
    setparameter(m, :ccm, :anomtable, transpose(anomtable))
    setparameter(m, :ccm, :temp, temp_forcing)

    run(m)
    return m[:ccm, :atmco2]
end



#Note on Parameter/Variable Values
#	deltat			=	timestep
#	Q10				=	Respiration Temperature sens.
#	Beta			=	Carbon Fertilization param.
#	Eta				=	Thermocline transfer velocity
#	atmco20			=	Initial concentration of CO2 (ppm)
#	temp_forcing	=	forcing of surface temperature (K)
#	co2_forcing 	=	forcing of co2 emissions
#
#	format for run_fortran_ccm = run_fortran_ccm(Climate_sensitivity, Soil_respiration, Carbon_fertilization, Thermocline_diffusion, Initial_CO2, Temp_Forcing, CO2_emissions_forcing)
#	Note: Climate Sensitivity parameter for run_fortran_ccm does not affect results (and is not included in mimi version)

#Run Mimi Verseion of ccm
m_atmco2 = mimiccm(1., 1.311, 0.502, 280., temp_forcing, co2_forcing);

#Run fortran version of ccm
init_fortran_ccm()
f_atmco2 = run_fortran_ccm(2., 1., 1.311, 0.502, 280., temp_forcing, co2_forcing);
fin_fortran_ccm()


#Test Precision
Precision = 1.0e-12

@test_approx_eq_eps maxabs(m_atmco2 .- f_atmco2) 0. Precision