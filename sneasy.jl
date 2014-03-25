using IAMF
include("doeclim.jl")
include("ccm.jl")
include("radforc.jl")

function getsneasy(nsteps=566)
	# ---------------------------------------------
	# Read data
	# ---------------------------------------------

	f_anomtable = readdlm("anomtable.txt")
	f_emissions = readdlm("emis_data_sep09.txt")
	f_nonco2forcing = readdlm("non_CO2_forcing.txt")

	# Timesteps
	deltat = 1.0

	# ---------------------------------------------
	# Create components
	# ---------------------------------------------

	c_doeclim = doeclimcomponent.doeclim(nsteps)
	c_ccm = ccmcomponent.ccm(nsteps)
	c_radforc = radforccomponent.radforc(nsteps)

	# ---------------------------------------------
	# Set parameters
	# ---------------------------------------------

	c_doeclim.p.t2co = 2.0
	c_doeclim.p.kappa = 1.1
	c_doeclim.p.deltat = deltat
	
	c_ccm.p.deltat = deltat
	c_ccm.p.Clim_sens=2.0
	c_ccm.p.Q10=1.311
	c_ccm.p.Beta=0.502
	c_ccm.p.Eta=17.722
	c_ccm.p.CO2_emissions=vec(f_emissions[:,2])
	c_ccm.p.anomtable=zeros(100,16000)
	for i=1:16000
	    c_ccm.p.anomtable[:,i] = f_anomtable[i,:]
	end
	
	c_radforc.p.other_forcing=vec(f_nonco2forcing[:,2])
	c_radforc.p.deltat=deltat

	# ---------------------------------------------
	# Connect parameters to variables
	# ---------------------------------------------

	c_doeclim.p.forcing = c_radforc.v.rf
	c_ccm.p.temp = c_doeclim.v.temp
	c_radforc.p.atmco2 = c_ccm.v.atmco2

	# ---------------------------------------------
	# Return model
	# ---------------------------------------------

	comps = [c_radforc,c_doeclim,c_ccm]
	return comps
end