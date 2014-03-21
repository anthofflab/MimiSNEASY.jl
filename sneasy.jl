include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")

# Read data
f_anomtable = readdlm("anomtable.txt")
f_emissions = readdlm("emis_data_sep09.txt")
f_nonco2forcing = readdlm("non_CO2_forcing.txt")

# Timesteps
nsteps = 200

pardoeclim = doeclim.doeclimpar(nsteps,1,3,0.5,zeros(nsteps))

parccm = ccm.ccmpar(nsteps,1,3,1.126598,0.2916273,172.2809,zeros(nsteps),vec(f_emissions[:,2]),zeros(100,16000))
parccm.CO2_emissions[:] = 7
for i=1:16000
    parccm.anomtable[:,i] = f_anomtable[i,:]
end

parradforc = radforc.radforcpar(nsteps,zeros(nsteps),vec(f_nonco2forcing[:,2]))

vardoeclim = doeclim.doeclimvar(pardoeclim)
varccm = ccm.ccmvar(parccm)
varradforc = radforc.radforcvar(parradforc)

# Connect parameters to vars
pardoeclim.forcing = varradforc.rf
parccm.temp = vardoeclim.temp
parradforc.atmco2 = varccm.atmco2

function runsneasy(pardoeclim,vardoeclim,parccm,varccm,parradforc,varradforc)

	doeclim.init(pardoeclim, vardoeclim)
	ccm.init(parccm, varccm)

	for t=1:nsteps
		radforc.timestep(parradforc, varradforc, t)
		doeclim.timestep(pardoeclim, vardoeclim, t)
		ccm.timestep(parccm, varccm, t)
	end
end
