using IAMF
include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")

function getsneasy(nsteps=566)
    m = Model()

    setindex(m, :time, nsteps)

    # ---------------------------------------------
    # Create components
    # ---------------------------------------------

    addcomponent(m, radforccomponent.radforc)
    addcomponent(m, doeclimcomponent.doeclim)
    addcomponent(m, ccmcomponent.ccm)    

    # ---------------------------------------------
    # Read data
    # ---------------------------------------------

    f_anomtable = readdlm("../data/anomtable.txt")
    f_emissions = readdlm("../data/emis_data_sep09.txt")
    f_nonco2forcing = readdlm("../data/non_CO2_forcing.txt")

    # Timesteps
    deltat = 1.0
    anomtable = zeros(100, 16000)
    for i=1:16000
        anomtable[:,i] = f_anomtable[i,:]
    end

    # ---------------------------------------------
    # Set parameters
    # ---------------------------------------------

    setparameter(m, :doeclim, :t2co, 2.0)
    setparameter(m, :doeclim, :kappa, 1.1)
    setparameter(m, :doeclim, :deltat, deltat)
    
    setparameter(m, :ccm, :deltat, deltat)
    setparameter(m, :ccm, :Clim_sens, 2.0)
    setparameter(m, :ccm, :Q10, 1.311)
    setparameter(m, :ccm, :Beta, 0.502)
    setparameter(m, :ccm, :Eta, 17.722)
    setparameter(m, :ccm, :CO2_emissions, vec(f_emissions[:,2]))
    setparameter(m, :ccm, :anomtable, anomtable)
    
    setparameter(m, :radforc, :other_forcing, vec(f_nonco2forcing[:,2]))
    setparameter(m, :radforc, :deltat, deltat)

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    bindparameter(m, :doeclim, :forcing, :radforc, :rf)
    bindparameter(m, :ccm, :temp, :doeclim)
    bindparameter(m, :radforc, :atmco2, :ccm)

    return m
end