include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")

# Read data
f_anomtable = readdlm("anomtable.txt")

# Timesteps
nsteps = 200

pardoeclim = doeclim.doeclimpar(nsteps,1,3,0.5,None)

parccm = ccm.ccmpar(nsteps,1,3,1.126598,0.2916273,172.2809,None,zeros(nsteps),zeros(100,16000))
parccm.CO2_emissions[:] = 7
for i=1:16000
    parccm.anomtable[:,i] = f_anomtable[i,:]
end

parradforc = radforc.radforcpar(nsteps,None,None)

vardoeclim = doeclim.doeclimvar(pardoeclim)
varccm = ccm.ccmvar(parccm)
varradforc = radforc.radforcvar(parradforc)

# Connect parameters to vars
pardoeclim.forcing = varradforc.rf
parccm.temp = vardoeclim.temp

doeclim.init(pardoeclim, vardoeclim)
ccm.init(parccm, varccm)

for t=1:nsteps
	radforc.timestep(parradforc, varradforc, t)
	doeclim.timestep(pardoeclim, vardoeclim, t)
	ccm.timestep(parccm, varccm, t)
end
