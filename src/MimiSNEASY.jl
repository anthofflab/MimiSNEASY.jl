module MimiSNEASY

using Mimi
using DataFrames

include("components/doeclim.jl")
include("components/ccm.jl")
include("components/radiativeforcing.jl")
include("components/rfco2.jl")

function getsneasy(;nsteps=736)
    m = Model()

    setindex(m, :time, nsteps)

    # ---------------------------------------------
    # Create components
    # ---------------------------------------------
    addcomponent(m, rfco2component.rfco2, :rfco2)
    addcomponent(m, radiativeforcingcomponent.radiativeforcing, :radiativeforcing)
    addcomponent(m, doeclimcomponent.doeclim, :doeclim)
    addcomponent(m, ccmcomponent.ccm, :ccm)

    # ---------------------------------------------
    # Read data
    # ---------------------------------------------

    f_anomtable = readdlm(joinpath(dirname(@__FILE__), "..", "data", "anomtable.txt"));
    rf_data = readtable(joinpath(dirname(@__FILE__), "..", "calibration", "data", "forcing_rcp85.txt"), separator = ' ', header=true);
    df = readtable(joinpath(dirname(@__FILE__), "..", "calibration", "data", "RCP85_EMISSIONS.csv"))
    rename!(df, :YEARS, :year);
    df = join(df, rf_data, on=:year, kind=:outer)
    df = DataFrame(year=df[:year], co2=df[:FossilCO2]+df[:OtherCO2], rf_aerosol=df[:aerosol_direct]+df[:aerosol_indirect], rf_other=df[:ghg_nonco2]+df[:volcanic]+df[:solar]+df[:other]);
    f_co2emissions = convert(Array, df[:co2]);
    f_rfaerosol = convert(Array, df[:rf_aerosol]);
    f_rfother = convert(Array, df[:rf_other]);

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

    setparameter(m, :radiativeforcing, :rf_aerosol, f_rfaerosol)
    #setparameter(m, :radiativeforcing, :rf_ch4, zeros(nsteps))
    setparameter(m, :radiativeforcing, :rf_other, f_rfother)
    setparameter(m, :radiativeforcing, :alpha, 1.)
    setparameter(m, :radiativeforcing, :deltat, deltat)

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    connectparameter(m, :doeclim, :forcing, :radiativeforcing, :rf)
    connectparameter(m, :ccm, :temp, :doeclim, :temp)
    connectparameter(m, :rfco2, :atmco2, :ccm, :atmco2)
    connectparameter(m, :radiativeforcing, :rf_co2, :rfco2, :rf_co2)

    return m
end

end
