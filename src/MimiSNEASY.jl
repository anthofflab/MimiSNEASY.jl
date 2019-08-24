module MimiSNEASY

using Mimi
using DataFrames
using DelimitedFiles
using SpecialFunctions

include("components/doeclim.jl")
include("components/ccm.jl")
include("components/radiativeforcing.jl")
include("components/rfco2.jl")

function getsneasy(;nsteps=736)
    m = Model()

    set_dimension!(m, :time, nsteps)
    set_dimension!(m, :anom_row, 100)
    set_dimension!(m, :anom_col, 16000)

    # ---------------------------------------------
    # Create components
    # ---------------------------------------------
    add_comp!(m, rfco2)
    add_comp!(m, radiativeforcing)
    add_comp!(m, doeclim)
    add_comp!(m, ccm)

    # ---------------------------------------------
    # Read data
    # ---------------------------------------------

    f_anomtable = readdlm(joinpath(dirname(@__FILE__), "..", "data", "anomtable.txt"));
    rf_data = readtable(joinpath(dirname(@__FILE__), "..", "calibration", "data", "forcing_rcp85.txt"), separator = ' ', header=true);
    df = readtable(joinpath(dirname(@__FILE__), "..", "calibration", "data", "RCP85_EMISSIONS.csv"))
    rename!(df, :YEARS => :year);
    df = join(df, rf_data, on=:year, kind=:outer)
    df = DataFrame(year=df.year, co2=df.FossilCO2+df.OtherCO2, rf_aerosol=df.aerosol_direct+df.aerosol_indirect, rf_other=df.ghg_nonco2+df.volcanic+df.solar+df.other);
    f_co2emissions = convert(Array, df.co2);
    f_rfaerosol = convert(Array, df.rf_aerosol);
    f_rfother = convert(Array, df.rf_other);

    # Timesteps
    deltat = 1.0
    anomtable = transpose(f_anomtable)

    # ---------------------------------------------
    # Set parameters
    # ---------------------------------------------

    set_param!(m, :doeclim, :t2co, 2.0)
    set_param!(m, :doeclim, :kappa, 1.1)
    set_param!(m, :doeclim, :deltat, deltat)

    set_param!(m, :ccm, :deltat, deltat)
    set_param!(m, :ccm, :Q10, 1.311)
    set_param!(m, :ccm, :Beta, 0.502)
    set_param!(m, :ccm, :Eta, 17.7)
    set_param!(m, :ccm, :atmco20, 280.)
    set_param!(m, :ccm, :CO2_emissions, f_co2emissions)
    set_param!(m, :ccm, :anomtable, anomtable)

    set_param!(m, :radiativeforcing, :rf_aerosol, f_rfaerosol)
    #set_param!(m, :radiativeforcing, :rf_ch4, zeros(nsteps))
    set_param!(m, :radiativeforcing, :rf_other, f_rfother)
    set_param!(m, :radiativeforcing, :alpha, 1.)
    set_param!(m, :radiativeforcing, :deltat, deltat)

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    connect_param!(m, :doeclim, :forcing, :radiativeforcing, :rf)
    connect_param!(m, :ccm, :temp, :doeclim, :temp)
    connect_param!(m, :rfco2, :CO₂, :ccm, :atmco2)
    connect_param!(m, :radiativeforcing, :rf_co2, :rfco2, :rf_co2)

    return m
end

end