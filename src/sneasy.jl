
using Mimi
using DataFrames

include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")

function getsneasy(;nsteps=736)
    m = Model()

    setindex(m, :time, nsteps)

    # ---------------------------------------------
    # Create components
    # ---------------------------------------------

    addcomponent(m, radforccomponent.radforc, :radforc)
    addcomponent(m, doeclimcomponent.doeclim, :doeclim)
    addcomponent(m, ccmcomponent.ccm, :ccm)

    # ---------------------------------------------
    # Read data
    # ---------------------------------------------

    f_anomtable = readdlm("../../sneasy/sneasy/anomtable.txt");
    f_nonco2forcing = readtable("../../sneasy/data/forcing_rcp85.txt", separator = ' ', header=true);
    df = readtable("../../sneasy/data/RCP85_EMISSIONS.csv");
    rename!(df, :YEARS, :year);
    df = DataFrame(year=df[:year], co2=df[:FossilCO2]+df[:OtherCO2]);
    f_co2emissions = convert(Array, df[:co2]);

    # Timesteps
    deltat = 1.0
    anomtable = transpose(f_anomtable)

    # ---------------------------------------------
    # Set parameters
    # ---------------------------------------------

    setparameter(m, :doeclim, :t2co, 2.0)
    setparameter(m, :doeclim, :kappa, 1.1)
    setparameter(m, :doeclim, :deltat, deltat)

    setparameter(m, :ccm, :deltat, deltat)
    setparameter(m, :ccm, :Q10, 1.311)
    setparameter(m, :ccm, :Beta, 0.502)
    setparameter(m, :ccm, :Eta, 17.7)
    setparameter(m, :ccm, :atmco20, 280.)
    setparameter(m, :ccm, :CO2_emissions, f_co2emissions)
    setparameter(m, :ccm, :anomtable, anomtable)

    setparameter(m, :radforc, :forcing_other, Vector(f_nonco2forcing[:other]))
    setparameter(m, :radforc, :forcing_volcanic, Vector(f_nonco2forcing[:volcanic]))
    setparameter(m, :radforc, :forcing_solar, Vector(f_nonco2forcing[:solar]))
    setparameter(m, :radforc, :forcing_ghg_nonco2, Vector(f_nonco2forcing[:ghg_nonco2]))
    setparameter(m, :radforc, :forcing_aerosol_direct, Vector(f_nonco2forcing[:aerosol_direct]))
    setparameter(m, :radforc, :forcing_aerosol_indirect, Vector(f_nonco2forcing[:aerosol_indirect]))
    setparameter(m, :radforc, :alpha, 1.)
    setparameter(m, :radforc, :deltat, deltat)

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    connectparameter(m, :doeclim, :forcing, :radforc, :rf)
    connectparameter(m, :ccm, :temp, :doeclim, :temp)
    connectparameter(m, :radforc, :atmco2, :ccm, :atmco2)

    return m
end