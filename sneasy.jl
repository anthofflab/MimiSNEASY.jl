using IAMF
include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")

function getsneasy(nsteps=566)
    # ---------------------------------------------
    # Read data
    # ---------------------------------------------

    f_anomtable = readdlm("data/anomtable.txt")
    f_emissions = readdlm("data/emis_data_sep09.txt")
    f_nonco2forcing = readdlm("data/non_CO2_forcing.txt")

    # Timesteps
    deltat = 1.0

    # ---------------------------------------------
    # Create components
    # ---------------------------------------------

    c_doeclim = doeclimcomponent.doeclim({:time=>nsteps})
    c_ccm = ccmcomponent.ccm({:time=>nsteps})
    c_radforc = radforccomponent.radforc({:time=>nsteps})

    # ---------------------------------------------
    # Set parameters
    # ---------------------------------------------

    c_doeclim.Parameters.t2co = 2.0
    c_doeclim.Parameters.kappa = 1.1
    c_doeclim.Parameters.deltat = deltat
    
    c_ccm.Parameters.deltat = deltat
    c_ccm.Parameters.Clim_sens=2.0
    c_ccm.Parameters.Q10=1.311
    c_ccm.Parameters.Beta=0.502
    c_ccm.Parameters.Eta=17.722
    c_ccm.Parameters.CO2_emissions=vec(f_emissions[:,2])
    c_ccm.Parameters.anomtable=zeros(100,16000)
    for i=1:16000
        c_ccm.Parameters.anomtable[:,i] = f_anomtable[i,:]
    end
    
    c_radforc.Parameters.other_forcing=vec(f_nonco2forcing[:,2])
    c_radforc.Parameters.deltat=deltat

    # ---------------------------------------------
    # Connect parameters to variables
    # ---------------------------------------------

    c_doeclim.Parameters.forcing = c_radforc.Variables.rf
    c_ccm.Parameters.temp = c_doeclim.Variables.temp
    c_radforc.Parameters.atmco2 = c_ccm.Variables.atmco2

    # ---------------------------------------------
    # Return model
    # ---------------------------------------------

    comps::Vector{ComponentState} = [c_radforc,c_doeclim,c_ccm]
    return comps
end