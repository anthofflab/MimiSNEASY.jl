module MimiSNEASY

using Mimi
using DataFrames
using DelimitedFiles
using SpecialFunctions

include("components/doeclim.jl")
include("components/ccm.jl")
include("components/radiativeforcing.jl")
include("components/rfco2.jl")

function getsneasy(;start_year::Int=1765, end_year::Int=2500)
    m = Model()

    set_dimension!(m, :time, start_year:end_year)
    set_dimension!(m, :anom_row, 100)
    set_dimension!(m, :anom_col, 16000)

    # Set number of model time steps.
    nsteps = length(start_year:end_year)

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

    # Read in RCP scenario and other data needed to run SNEASY.
    df = readtable(joinpath(dirname(@__FILE__), "..", "calibration", "data", "RCP85_EMISSIONS.csv"))
    rf_data = readtable(joinpath(dirname(@__FILE__), "..", "calibration", "data", "forcing_rcp85.txt"), separator = ' ', header=true)
    f_anomtable = readdlm(joinpath(dirname(@__FILE__), "..", "data", "anomtable.txt"));

    # Get RCP year indices based on user-specified time horizon to run model.
    start_index, end_index = findall((in)([start_year, end_year]), collect(1765:2500))

    # Clean up model data.
    rename!(df, :YEARS => :year);
    df = join(df, rf_data, on=:year, kind=:outer)
    df = DataFrame(year=df.year, co2=df.FossilCO2+df.OtherCO2, rf_aerosol=df.aerosol_direct+df.aerosol_indirect, rf_other=df.ghg_nonco2+df.volcanic+df.solar+df.other);

    # Get specific input variables indexed to user-specified model time horizon.
    f_co2emissions = convert(Array, df.co2)[start_index:end_index]
    f_rfaerosol = convert(Array, df.rf_aerosol)[start_index:end_index]
    f_rfother = convert(Array, df.rf_other)[start_index:end_index]

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
    # TODO Find better values for the following two parameters
    set_param!(m, :ccm, :oxidised_CH₄_to_CO₂, zeros(nsteps))
    set_param!(m, :rfco2, :N₂O, fill(272.95961, nsteps))

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
