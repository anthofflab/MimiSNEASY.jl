#  DOECLIM:  Diffusion Ocean Energy balance CLIMate model
#------------------------------------------------------------------------------
# Simple climate model DOECLIM
#
# calculates sea surface and land air temperature response to radiative forcing
# based on an energy balance model with 1-D diffusion ocean
#
# Constructed by Elmar Kriegler (EK),
# Potsdam Institute for Climate Impact Research
# Date: 06.02.2005
#
# References for the historical forcing values can be found in (Reference EK05):
# Kriegler, E (2005) Imprecise probability analysis for integrated assessment
# of climate change. Ph.D. thesis. University of Potsdam, 256 pp.
# opus.kobv.de/ubp/volltexte/2005/561/
#
# Model equations are described in EK05 and (Reference TK07):
# Tanaka, K, Kriegler, E, Bruckner, T, Hooss, G, Knorr, W, Raddatz, T (2007)
# Aggregated carbon cycle, atmospheric chemistry, and climate model (ACC2):
# Description of the forward and inverse modes, Reports on Earth System Science ! 40/2007,
# Max Planck Institute for Meteorology, Hamburg, 199 pp.
# www.mpimet.mpg.de/fileadmin/publikationen/Reports/BzE_40.pdf
#
#==============================================================================
#
# Updates:
# 22.05.2007 Hammer-Hollingsworth numerical correction included (EK)
# 23.05.2007 Ocean heat uptake added (EK)
# 12.02.2008 Translated to Fortran90 (Marlos Goes <mpg14@psu.edu>)
# 15.08.2009 Written as Fortran90 module (Brian Tuttle <btuttle@psu.edu>)
#  
#==============================================================================
#
# Global Parameters:
#   ak      slope coeff. for land-sea heat exchange
#   bk      inters. coeff. for land-sea heat exch.
#   bsi     marine air warming enhancement
#   cal     heat cap. of land-troposph. system
#   cas     heat cap. of ocean ML-troposph.
#   csw     specific heat capacity of 1m^3 seawater [Wa/m^3/K]
#   deltat  time step size [years]
#   flnd    land fraction
#   fso     ocean frac. area below 60m
#   kcon    conversion factor [cm2/s->m2/a]
#   q2co    2xCo2 forcing increase [W/m^2]
#   rlam    clim sens. over land enhancement
#   zbot    depth of interior ocean
#   
#   temp_landair:       land air temperature anomaly (K)
#   temp_sst:           sea surface temperature anomaly (K)
#   heat_mixed:         mixed layer heat anomaly (10^22 J)
#   heat_interior:      interior ocean heat anomaly (10^22 J)
#   heatflux_mixed:     heat uptake of the mixed layer (W/m^2)
#   heatflux_interior:  heat uptake of the interior ocean (W/m^2)
#
#==============================================================================
module doeclim

const ak   = 0.31
const bk   = 1.59
const bsi  = 1.3
const cal  = 0.52
const cas  = 7.80
const csw  = 0.13
const flnd = 0.29
const fso  = 0.95
const kcon = 3155.8
const q2co = 3.7
const rlam = 1.43
const zbot = 4000
const earth_area = 5100656e8      # [m^2]
const secs_per_Year = 31556926.0

type doeclimpar
	nsteps::Int
	deltat::Float64
	# climate sensitivity to 2xCO2 (K); default = 3
	t2co::Float64              
	# vertical ocean diffusivity (cm^2 s^-1); default = 0.55
	kappa::Float64             
	forcing::Vector{Float64}
end

type doeclimvar
	temp::Vector{Float64}
	temp_landair::Vector{Float64}
	temp_sst::Vector{Float64}
	heat_mixed::Vector{Float64}
	heat_interior::Vector{Float64}
	heatflux_mixed::Vector{Float64}
	heatflux_interior::Vector{Float64}
	Ker::Vector{Float64}
	IB::Matrix{Float64}
	Adoe::Matrix{Float64}
	taucfl::Float64
	taukls::Float64
	taucfs::Float64
	tauksl::Float64
	taudif::Float64
	taubot::Float64
	powtoheat::Float64

	KT0::Vector{Float64}
	KTA1::Vector{Float64}
	KTB1::Vector{Float64}
	KTA2::Vector{Float64}
	KTB2::Vector{Float64}
	KTA3::Vector{Float64}
	KTB3::Vector{Float64}

	Cdoe::Matrix{Float64}
	Baux::Matrix{Float64}

	function doeclimvar(p::doeclimpar)
		vars = new(
			zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(2,2),
      		zeros(2,2),
      		0,
      		0,
      		0,
      		0,
      		0,
      		0,
      		0,
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(p.nsteps),
      		zeros(2,2),
      		zeros(2,2))
    	return vars
  	end
end

function init(p::doeclimpar, v::doeclimvar)
	# DEPENDENT MODEL PARAMETERS
	ocean_area = (1.0-flnd)*earth_area

	v.powtoheat = ocean_area*secs_per_Year / 1e22

	cnum = rlam*flnd + bsi * (1.0-flnd)

	cden = rlam * flnd - ak *(rlam-bsi)

	# vertical diffusivity in [m^2/a]

    keff = kcon * p.kappa

	# climate feedback strength over land

    cfl = flnd *cnum/cden*q2co/p.t2co-bk*(rlam-bsi)/cden

	# climate feedback strength over ocean

    cfs = (rlam * flnd - ak / (1.-flnd) * (rlam-bsi)) * cnum / cden * q2co / p.t2co + rlam * flnd / (1.-flnd) * bk * (rlam - bsi) / cden

	# land-sea heat exchange coefficient

    kls = bk * rlam * flnd / cden - ak * flnd * cnum / cden * q2co / p.t2co

	# interior ocean warming time scale

    v.taubot = zbot^2 / keff

	# ocean heat diff. time scale

    v.taudif = cas^2 / csw^2 * pi / keff

	# ocean response time scale

    v.taucfs = cas / cfs

	# land response time scale

    v.taucfl = cal / cfl

	# sea-land heat exchange time scale

    v.tauksl  = (1.-flnd) * cas / kls

	# land-sea heat exchange time scale

    v.taukls  = flnd * cal / kls

	# Zeroth Order

    v.KT0[p.nsteps] = 4-2*sqrt(2.)

	# First Order

    v.KTA1[p.nsteps] = -8*exp(-v.taubot/p.deltat) + 4*sqrt(2.)*exp(-0.5*v.taubot/p.deltat)

    v.KTB1[p.nsteps] = 4*sqrt(pi*v.taubot/p.deltat) * (1+erf(sqrt(0.5*v.taubot/p.deltat)) - 2*erf(sqrt(v.taubot/p.deltat)))

	# Second order

    v.KTA2[p.nsteps] =  8*exp(-4.*v.taubot/p.deltat) - 4*sqrt(2.)*exp(-2.*v.taubot/p.deltat)

	v.KTB2[p.nsteps] = -8*sqrt(pi*v.taubot/p.deltat) * (1.+ erf(sqrt(2.*v.taubot/p.deltat)) - 2.*erf(2.*sqrt(v.taubot/p.deltat)) )

	# Third Order

    v.KTA3[p.nsteps] = -8.*exp(-9.*v.taubot/p.deltat) + 4*sqrt(2.)*exp(-4.5*v.taubot/p.deltat)

    v.KTB3[p.nsteps] = 12.*sqrt(pi*v.taubot/p.deltat) * (1 +erf(sqrt(4.5*v.taubot/p.deltat)) - 2.*erf(3.*sqrt(v.taubot/p.deltat)) )

	# Hammer and Hollingsworth correction (Equation 2.3.27, TK07):
	# Switched on (To switch off, comment out lines below)
    v.Cdoe[1,1] = 1./v.taucfl^2+1./v.taukls^2+2./v.taucfl/v.taukls+bsi/v.taukls/v.tauksl*(p.deltat^2/12.)
    v.Cdoe[1,2] = -bsi/v.taukls^2-bsi/v.taucfl/v.taukls-bsi/v.taucfs/v.taukls-bsi^2/v.taukls/v.tauksl*(p.deltat^2/12.)
    v.Cdoe[2,1] = -bsi/v.tauksl^2-1./v.taucfs/v.tauksl-1./v.taucfl/v.tauksl-1./v.taukls/v.tauksl*(p.deltat^2/12.)
    v.Cdoe[2,2] =  1./v.taucfs^2+bsi^2/v.tauksl^2+2.*bsi/v.taucfs/v.tauksl+bsi/v.taukls/v.tauksl*(p.deltat^2/12.)

	# Matrices of difference equation system B*T(i+1) = Q(i) + A*T(i)
	# T = (TL,TO)
	# (Equation A.27, EK05, or Equations 2.3.24 and 2.3.27, TK07)
    v.Baux[1,1] = 1. + p.deltat/(2.*v.taucfl) + p.deltat/(2.*v.taukls)+v.Cdoe[1,1]
    v.Baux[1,2] = -p.deltat/(2.*v.taukls)*bsi+v.Cdoe[1,2]
    v.Baux[2,1] = -p.deltat/(2.*v.tauksl)+v.Cdoe[2,1]
    v.Baux[2,2] = 1. + p.deltat/(2.*v.taucfs) + p.deltat/(2.*v.tauksl)*bsi + 2.*fso*sqrt(p.deltat/v.taudif)+v.Cdoe[2,2]

    v.IB=inv(v.Baux)

	for i = 1:p.nsteps-1

		# Zeroth Order

      	v.KT0[i] = 4*sqrt(float64(p.nsteps+1-i)) - 2.*sqrt(float64(p.nsteps+2-i)) - 2*sqrt(float64(p.nsteps-i))

		# First Order

		v.KTA1[i] = -8*sqrt(float64(p.nsteps+1-i)) * exp(-v.taubot/p.deltat/(p.nsteps+1-i)) + 4*sqrt(float64(p.nsteps+2-i)) *exp(-v.taubot/p.deltat/(p.nsteps+2-i)) + 4*sqrt(float64(p.nsteps-i)) *exp(-v.taubot/p.deltat/(p.nsteps-i))

    	v.KTB1[i] = 4*sqrt(pi*v.taubot/p.deltat) * (erf(sqrt(v.taubot/p.deltat/(p.nsteps-i))) + erf(sqrt(v.taubot/p.deltat/(p.nsteps+2-i))) - 2*erf(sqrt(v.taubot/p.deltat/(p.nsteps+1-i))) )

		# Second Order

		v.KTA2[i] =  8.*sqrt(float64(p.nsteps+1-i)) * exp(-4.*v.taubot/p.deltat/(p.nsteps+1-i))- 4.*sqrt(float64(p.nsteps+2-i))*exp(-4.*v.taubot/p.deltat/(p.nsteps+2-i))- 4.*sqrt(float64(p.nsteps-i)) * exp(-4.*v.taubot/p.deltat/(p.nsteps-i))

		v.KTB2[i] = -8.*sqrt(pi*v.taubot/p.deltat) * (erf(2.*sqrt(v.taubot/p.deltat/(float64(p.nsteps-i)))) + erf(2.*sqrt(v.taubot/p.deltat/float64(p.nsteps+2-i))) -       2.*erf(2.*sqrt(v.taubot/p.deltat/float64(p.nsteps+1-i))) )

		# Third Order

		v.KTA3[i] = -8.*sqrt(float64(p.nsteps+1-i)) *exp(-9.*v.taubot/p.deltat/(p.nsteps+1.-i)) + 4.*sqrt(float64(p.nsteps+2-i))*exp(-9.*v.taubot/p.deltat/(p.nsteps+2.-i)) + 4.*sqrt(float64(p.nsteps-i))*exp(-9.*v.taubot/p.deltat/(p.nsteps-i))

		v.KTB3[i] = 12.*sqrt(pi*v.taubot/p.deltat) * (erf(3.*sqrt(v.taubot/p.deltat/(p.nsteps-i))) + erf(3.*sqrt(v.taubot/p.deltat/(p.nsteps+2-i))) - 2.*erf(3.*sqrt(v.taubot/p.deltat/(p.nsteps+1-i))) )
	end

	for i=1:p.nsteps
		v.Ker[i] = v.KT0[i]+v.KTA1[i]+v.KTB1[i]+v.KTA2[i]+v.KTB2[i]+v.KTA3[i]+v.KTB3[i]
	end

	v.Adoe[1,1] = 1 - p.deltat/(2.*v.taucfl) - p.deltat/(2.*v.taukls) + v.Cdoe[1,1]
    v.Adoe[1,2] =  p.deltat/(2.*v.taukls)*bsi + v.Cdoe[1,2]
    v.Adoe[2,1] =  p.deltat/(2.*v.tauksl) + v.Cdoe[2,1]
    v.Adoe[2,2] = 1 - p.deltat/(2.*v.taucfs) - p.deltat/(2.*v.tauksl)*bsi + v.Ker[p.nsteps]*fso*sqrt(p.deltat/v.taudif) + v.Cdoe[2,2]
end

function timestep(p::doeclimpar, s::doeclimvar, n::Int)
#  ==========================================================================
# | Simple climate model DOECLIM
# |
# | calculates sea surface and land air temperature response to radiative 
# | forcing based on an energy balance model with 1-D diffusion ocean
# |
# | *** computes single time step ***
# | *** initialize with init_doeclim ***
# | *** then iterate this function ***
# |
# | Input:
# |       n:        current time step
# |       forcing:  global radiative forcing (top of atmosphere) (W/m^2)
# |
# | Output:
# |       temp: global mean temperature anomaly (K), relative to preindustrial
# |
# | Assumptions: 
# |       land surface temperature = land air temperature
# |       mixed layer temperature  = sea surface temperatures 
# |                                = marine air temperature divided by bsi
#  ==========================================================================#
	DTE1 = s.temp_landair
	DTE2 = s.temp_sst

	# assume land and ocean forcings are equal to global forcing
	QL = p.forcing
    Q0 = p.forcing

	if n>2
		DelQL = QL[n] - QL[n-1]
     	DelQ0 = Q0[n] - Q0[n-1]

		# Assumption: linear forcing change between n and n+1
		QC1 = (DelQL/cal*(1./s.taucfl+1./s.taukls)-bsi*DelQ0/cas/s.taukls)
		QC2 = (DelQ0/cas*(1./s.taucfs+bsi/s.tauksl)-DelQL/cal/s.tauksl)

		QC1 = QC1 * p.deltat^2/12.
		QC2 = QC2 * p.deltat^2/12.
		# -------------------------- INITIAL CONDITIONS ------------------------
		# Initialization of temperature and forcing vector:
		# Factor 1/2 in front of Q in Equation A.27, EK05, and Equation 2.3.27, TK07 is a typo!
		# Assumption: linear forcing change between n and n+1
		DQ1 = 0.5*p.deltat/cal*(QL[n]+QL[n-1])
		DQ2 = 0.5*p.deltat/cas*(Q0[n]+Q0[n-1])
		DQ1 = DQ1 + QC1
		DQ2 = DQ2 + QC2

		# -------------- SOLVE MODEL ------------------------------------
		# Calculate temperatures
		DPAST1 = 0.0
		DPAST2 = 0.0
		for i=1:n-1
			DPAST2 = DPAST2+DTE2[i]*s.Ker[p.nsteps-n+i]
		end

		DPAST2 = DPAST2*fso * sqrt(p.deltat/s.taudif)

		DTEAUX1 = s.Adoe[1,1]*DTE1[n-1]+s.Adoe[1,2]*DTE2[n-1]
		DTEAUX2 = s.Adoe[2,1]*DTE1[n-1]+s.Adoe[2,2]*DTE2[n-1]

		DTE1[n] = s.IB[1,1]*(DQ1+DPAST1+DTEAUX1)+s.IB[1,2]*(DQ2+DPAST2+DTEAUX2)
		DTE2[n] = s.IB[2,1]*(DQ1+DPAST1+DTEAUX1)+s.IB[2,2]*(DQ2+DPAST2+DTEAUX2)
	else
		# handle initial values

		DTE1[1] = 0.0
		DTE2[1] = 0.0

		DelQL = QL[2] - QL[1]
		DelQ0 = Q0[2] - Q0[1]

		QC1 = DelQL/cal*(1./s.taucfl+1./s.taukls)-bsi*DelQ0/cas/s.taukls
		QC2 = DelQ0/cas*(1./s.taucfs+bsi/s.tauksl)-DelQL/cal/s.tauksl
		QC1 = QC1*p.deltat^2/12.
		QC2 = QC2*p.deltat^2/12.

		DQ1 = 0.5*p.deltat/cal*(QL[2]+QL[1])
		DQ2 = 0.5*p.deltat/cas*(Q0[2]+Q0[1])
		DQ1= DQ1 + QC1
		DQ2= DQ2 + QC2

		DTEAUX1=s.Adoe[1,1]*DTE1[1]+s.Adoe[1,2]*DTE2[1]
		DTEAUX2=s.Adoe[2,1]*DTE1[1]+s.Adoe[2,2]*DTE2[1]

		DTE1[2] = s.IB[1,1]*(DQ1 +DTEAUX1)+s.IB[1,2]*(DQ2+DTEAUX2)
		DTE2[2] = s.IB[2,1]*(DQ1 +DTEAUX1)+s.IB[2,2]*(DQ2+DTEAUX2)
	end
	
	s.temp[n] = flnd*s.temp_landair[n] + (1.-flnd)*bsi*s.temp_sst[n]

	# Calculate ocean heat uptake [W/m^2]
	# heatflux(n) captures in the heat flux in the period between n-1 and n
	# Numerical implementation of Equation 2.7, EK05, or Equation 2.3.13, TK07)
	# ------------------------------------------------------------------------


	if n>1
		s.heatflux_mixed[n] = cas*( DTE2[n]-DTE2[n-1] )

		for i=1:n-1
        	s.heatflux_interior[n] = s.heatflux_interior[n]+DTE2[i]*s.Ker[p.nsteps-n+1+i]
		end
		s.heatflux_interior[n] = cas*fso/sqrt(s.taudif*p.deltat)*(2.*DTE2[n] - s.heatflux_interior[n])

		s.heat_mixed[n] = s.heat_mixed[n-1] +s.heatflux_mixed[n] *(s.powtoheat*p.deltat)

		s.heat_interior[n] = s.heat_interior[n-1] + s.heatflux_interior[n] * (fso*s.powtoheat*p.deltat)

	else
		# handle initial values
		s.heatflux_mixed[1] = 0.0
		s.heatflux_interior[1] = 0.0

		s.heat_mixed[1] = 0.0
		s.heat_interior[1] = 0.0
	end

end
end
