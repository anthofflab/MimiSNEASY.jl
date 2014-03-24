using Iamf
include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")
import radforc
import ccm
import doeclim

function getsneasy()
	# Read data
	f_anomtable = readdlm("anomtable.txt")
	f_emissions = readdlm("emis_data_sep09.txt")
	f_nonco2forcing = readdlm("non_CO2_forcing.txt")

	# Timesteps
	nsteps = 566
	deltat = 1.0

	# Create parameter objects for each component
	pardoeclim = doeclim.doeclimpar(nsteps,deltat)
	pardoeclim.t2co = 2.0
	pardoeclim.kappa = 1.1

	parccm = ccm.ccmpar(nsteps,deltat)
	parccm.Clim_sens=2.0
	parccm.Q10=1.311
	parccm.Beta=0.502
	parccm.Eta=17.722
	parccm.CO2_emissions=vec(f_emissions[:,2])
	parccm.anomtable=zeros(100,16000)
	for i=1:16000
	    parccm.anomtable[:,i] = f_anomtable[i,:]
	end

	parradforc = radforc.radforcpar(nsteps,deltat)
	parradforc.other_forcing=vec(f_nonco2forcing[:,2])

	# Create variable objects for each component
	vardoeclim = doeclim.doeclimvar(pardoeclim)
	varccm = ccm.ccmvar(parccm)
	varradforc = radforc.radforcvar(parradforc)

	# Connect parameters to vars
	pardoeclim.forcing = varradforc.rf
	parccm.temp = vardoeclim.temp
	parradforc.atmco2 = varccm.atmco2

	c1 = ComponentInstance(parradforc,varradforc)
	c2 = ComponentInstance(pardoeclim,vardoeclim)
	c3 = ComponentInstance(parccm,varccm)

	comps = [c1,c2,c3]
	return comps
end