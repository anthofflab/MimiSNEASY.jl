#!  DOECLIM:  Diffusion Ocean Energy balance CLIMate model
#!
#!  Copyright (C) 2007 E. Kriegler
#!
#!  This program is free software; you can redistribute it and/or modify
#!  it under the terms of the GNU General Public License as published by
#!  the Free Software Foundation; either version 2 of the License, or
#!  (at your option) any later version.
#!
#!  This program is distributed in the hope that it will be useful,
#!  but WITHOUT ANY WARRANTY; without even the implied warranty of
#!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#!  GNU General Public License for more details.
#!
#!  You should have received a copy of the GNU General Public License
#!  along with this program; if not, write to the Free Software
#!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#!
#!------------------------------------------------------------------------------
#! Simple climate model DOECLIM
#!
#! calculates sea surface and land air temperature response to radiative forcing
#! based on an energy balance model with 1-D diffusion ocean
#!
#! Constructed by Elmar Kriegler (EK),
#! Potsdam Institute for Climate Impact Research
#! Date: 06.02.2005
#!
#! References for the historical forcing values can be found in (Reference EK05):
#! Kriegler, E (2005) Imprecise probability analysis for integrated assessment
#! of climate change. Ph.D. thesis. University of Potsdam, 256 pp.
#! opus.kobv.de/ubp/volltexte/2005/561/
#!
#! Model equations are described in EK05 and (Reference TK07):
#! Tanaka, K, Kriegler, E, Bruckner, T, Hooss, G, Knorr, W, Raddatz, T (2007)
#! Aggregated carbon cycle, atmospheric chemistry, and climate model (ACC2):
#! Description of the forward and inverse modes, Reports on Earth System Science ! 40/2007,
#! Max Planck Institute for Meteorology, Hamburg, 199 pp.
#! www.mpimet.mpg.de/fileadmin/publikationen/Reports/BzE_40.pdf
#!
#!==============================================================================
#!
#! Updates:
#! 22.05.2007 Hammer-Hollingsworth numerical correction included (EK)
#! 23.05.2007 Ocean heat uptake added (EK)
#! 12.02.2008 Translated to Fortran90 (Marlos Goes <mpg14@psu.edu>)
#! 15.08.2009 Written as Fortran90 module (Brian Tuttle <btuttle@psu.edu>)
#!  
#!==============================================================================
#!
#! Global Parameters:
#!   ak      slope coeff. for land-sea heat exchange
#!   bk      inters. coeff. for land-sea heat exch.
#!   bsi     marine air warming enhancement
#!   cal     heat cap. of land-troposph. system
#!   cas     heat cap. of ocean ML-troposph.
#!   csw     specific heat capacity of 1m^3 seawater [Wa/m^3/K]
#!   deltat  time step size [years]
#!   flnd    land fraction
#!   fso     ocean frac. area below 60m
#!   kcon    conversion factor [cm2/s->m2/a]
#!   q2co    2xCo2 forcing increase [W/m^2]
#!   rlam    clim sens. over land enhancement
#!   zbot    depth of interior ocean
#!   
#!   temp_landair:       land air temperature anomaly (K)
#!   temp_sst:           sea surface temperature anomaly (K)
#!   heat_mixed:         mixed layer heat anomaly (10^22 J)
#!   heat_interior:      interior ocean heat anomaly (10^22 J)
#!   heatflux_mixed:     heat uptake of the mixed layer (W/m^2)
#!   heatflux_interior:  heat uptake of the interior ocean (W/m^2)
#!
#!==============================================================================
module doeclim

type doeclimvar
  nsteps
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

  function doeclimvar(timesteps)
    vars = new(
      timesteps,
      Array(Float64,timesteps),
      Array(Float64,timesteps),
      Array(Float64,timesteps),
      Array(Float64,timesteps),
      Array(Float64,timesteps),
      Array(Float64,timesteps),
      Array(Float64,timesteps),
      Array(Float64,2,2),
      Array(Float64,2,2),
      0,
      0,
      0,
      0,
      0,
      0,
      0)
    return vars
  end
end

function init_doeclim_arrays(doeclimvars)
#  ==========================================================================
# |  This routine allocates and initializes global arrays for DOECLIM.  It   |
# |  is separate from init_doeclim_parameters so that the parameters can be  |
# |  changed without reinitializing the arrays.                              |
#  ==========================================================================
  # Initialize global arrays to zero.
  doeclimvars.temp_landair[:]  = 0.0
  doeclimvars.temp_sst[:] = 0.0
  doeclimvars.heat_mixed[:] = 0.0
  doeclimvars.heat_interior[:] = 0.0
  doeclimvars.heatflux_mixed[:] = 0.0
  doeclimvars.heatflux_interior[:] = 0.0
end

function init_doeclim_parameters(doeclimvars, deltat, t2co, kappa)
#  =========================================================================
# |  Initialize variables for DOECLIM.                                      |
# |                                                                         |
# |  Input parameters:                                                      |
# |     t2co:   climate sensitivity to 2xCO2 (K); default = 3               |
# |     kappa:  vertical ocean diffusivity (cm^2 s^-1); default = 0.55      |
# |     nsteps: number of steps (length of forcing and response vectors)    |
# |                                                                         |
#  =========================================================================
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
	const secs_per_Year = 31556926

	nsteps = doeclimvars.nsteps
	KT0 = Array(Float64,nsteps)
	KTA1 = Array(Float64,nsteps)
	KTB1 = Array(Float64,nsteps)
	KTA2 = Array(Float64,nsteps)
	KTB2 = Array(Float64,nsteps)
	KTA3 = Array(Float64,nsteps)
	KTB3 = Array(Float64,nsteps)

	Cdoe = Array(Float64, 2, 2)
	Baux = Array(Float64, 2, 2)

#    real(DP) :: cden

	# DEPENDENT MODEL PARAMETERS
	ocean_area = (1.0-flnd)*earth_area

	doeclimvars.powtoheat = ocean_area*secs_per_Year / 1e22

	cnum = rlam*flnd + bsi * (1.0-flnd)

	cden = rlam * flnd - ak *(rlam-bsi)


	# vertical diffusivity in [m^2/a]

    keff = kcon * kappa

	# climate feedback strength over land

    cfl = flnd *cnum/cden*q2co/t2co-bk*(rlam-bsi)/cden

	# climate feedback strength over ocean

    cfs = (rlam * flnd - ak / (1.-flnd) * (rlam-bsi)) * cnum / cden * q2co / t2co + rlam * flnd / (1.-flnd) * bk * (rlam - bsi) / cden

	# land-sea heat exchange coefficient

    kls = bk * rlam * flnd / cden - ak * flnd * cnum / cden * q2co / t2co

	# interior ocean warming time scale

    doeclimvars.taubot = zbot^2 / keff

	# ocean heat diff. time scale

    doeclimvars.taudif = cas^2 / csw^2 * pi / keff

	# ocean response time scale

    doeclimvars.taucfs = cas / cfs

	# land response time scale

    doeclimvars.taucfl = cal / cfl

	# sea-land heat exchange time scale

    doeclimvars.tauksl  = (1.-flnd) * cas / kls

	# land-sea heat exchange time scale

    doeclimvars.taukls  = flnd * cal / kls


	# Zeroth Order

    KT0[nsteps] = 4-2*sqrt(2.)

	# First Order

    KTA1[nsteps] = -8*exp(-doeclimvars.taubot/deltat) + 4*sqrt(2.)*exp(-0.5*doeclimvars.taubot/deltat)

    KTB1[nsteps] = 4*sqrt(pi*doeclimvars.taubot/deltat) * (1+erf(sqrt(0.5*doeclimvars.taubot/deltat)) - 2*erf(sqrt(doeclimvars.taubot/deltat)))

	# Second order

    KTA2[nsteps] =  8*exp(-4.*doeclimvars.taubot/deltat) - 4*sqrt(2.)*exp(-2.*doeclimvars.taubot/deltat)

	KTB2[nsteps] = -8*sqrt(pi*doeclimvars.taubot/deltat) * (1.+ erf(sqrt(2.*doeclimvars.taubot/deltat)) - 2.*erf(2.*sqrt(doeclimvars.taubot/deltat)) )

	# Third Order

    KTA3[nsteps] = -8.*exp(-9.*doeclimvars.taubot/deltat) + 4*sqrt(2.)*exp(-4.5*doeclimvars.taubot/deltat)

    KTB3[nsteps] = 12.*sqrt(pi*doeclimvars.taubot/deltat) * (1 +erf(sqrt(4.5*doeclimvars.taubot/deltat)) - 2.*erf(3.*sqrt(doeclimvars.taubot/deltat)) )

	# %Hammer and Hollingsworth correction (Equation 2.3.27, TK07):
	# %Switched on (To switch off, comment out lines below)
    Cdoe[1,1] = 1./doeclimvars.taucfl^2+1./doeclimvars.taukls^2+2./doeclimvars.taucfl/doeclimvars.taukls+bsi/doeclimvars.taukls/doeclimvars.tauksl
    Cdoe[1,2] = -bsi/doeclimvars.taukls^2-bsi/doeclimvars.taucfl/doeclimvars.taukls-bsi/doeclimvars.taucfs/doeclimvars.taukls-bsi^2/doeclimvars.taukls/doeclimvars.tauksl
    Cdoe[2,1] = -bsi/doeclimvars.tauksl^2-1./doeclimvars.taucfs/doeclimvars.tauksl-1./doeclimvars.taucfl/doeclimvars.tauksl-1./doeclimvars.taukls/doeclimvars.tauksl
    Cdoe[2,2] =  1./doeclimvars.taucfs^2+bsi^2/doeclimvars.tauksl^2+2.*bsi/doeclimvars.taucfs/doeclimvars.tauksl+bsi/doeclimvars.taukls/doeclimvars.tauksl
    Cdoe = Cdoe*(deltat^2/12.)

	#%------------------------------------------------------------------
	#% Matrices of difference equation system B*T(i+1) = Q(i) + A*T(i)
	#% T = (TL,TO)
	#% (Equation A.27, EK05, or Equations 2.3.24 and 2.3.27, TK07)
    Baux[1,1] = 1. + deltat/(2.*doeclimvars.taucfl) + deltat/(2.*doeclimvars.taukls)
    Baux[1,2] = -deltat/(2.*doeclimvars.taukls)*bsi
    Baux[2,1] = -deltat/(2.*doeclimvars.tauksl)
    Baux[2,2] = 1. + deltat/(2.*doeclimvars.taucfs) + deltat/(2.*doeclimvars.tauksl)*bsi + 2.*fso*sqrt(deltat/doeclimvars.taudif)
    Baux = Baux+Cdoe

	# Calculate inverse of B
    #call migs(Baux,2,IBaux)
    #IB[:,:]=IBaux(:,:)
    doeclimvars.IB[:,:]=inv(Baux)[:,:]

	for i = 1:nsteps-1

		# Zeroth Order

      	KT0[i] = 4*sqrt(float64(nsteps+1-i)) - 2.*sqrt(float64(nsteps+2-i)) - 2*sqrt(float64(nsteps-i))

		# First Order

		KTA1[i] = -8*sqrt(float64(nsteps+1-i)) * exp(-doeclimvars.taubot/deltat/(nsteps+1-i)) + 4*sqrt(float64(nsteps+2-i)) *exp(-doeclimvars.taubot/deltat/(nsteps+2-i)) + 4*sqrt(float64(nsteps-i)) *exp(-doeclimvars.taubot/deltat/(nsteps-i))

    	KTB1[i] = 4*sqrt(pi*doeclimvars.taubot/deltat) * (erf(sqrt(doeclimvars.taubot/deltat/(nsteps-i))) + erf(sqrt(doeclimvars.taubot/deltat/(nsteps+2-i))) - 2*erf(sqrt(doeclimvars.taubot/deltat/(nsteps+1-i))) )

		# Second Order

		KTA2[i] =  8.*sqrt(float64(nsteps+1-i)) * exp(-4.*doeclimvars.taubot/deltat/(nsteps+1-i))- 4.*sqrt(float64(nsteps+2-i))*exp(-4.*doeclimvars.taubot/deltat/(nsteps+2-i))- 4.*sqrt(float64(nsteps-i)) * exp(-4.*doeclimvars.taubot/deltat/(nsteps-i))

		KTB2[i] = -8.*sqrt(pi*doeclimvars.taubot/deltat) * (erf(2.*sqrt(doeclimvars.taubot/deltat/(float64(nsteps-i)))) + erf(2.*sqrt(doeclimvars.taubot/deltat/float64(nsteps+2-i))) -       2.*erf(2.*sqrt(doeclimvars.taubot/deltat/float64(nsteps+1-i))) )

		# Third Order

		KTA3[i] = -8.*sqrt(float64(nsteps+1-i)) *exp(-9.*doeclimvars.taubot/deltat/(nsteps+1.-i)) + 4.*sqrt(float64(nsteps+2-i))*exp(-9.*doeclimvars.taubot/deltat/(nsteps+2.-i)) + 4.*sqrt(float64(nsteps-i))*exp(-9.*doeclimvars.taubot/deltat/(nsteps-i))

		KTB3[i] = 12.*sqrt(pi*doeclimvars.taubot/deltat) * (erf(3.*sqrt(doeclimvars.taubot/deltat/(nsteps-i))) + erf(3.*sqrt(doeclimvars.taubot/deltat/(nsteps+2-i))) - 2.*erf(3.*sqrt(doeclimvars.taubot/deltat/(nsteps+1-i))) )
	end

	doeclimvars.Ker = KT0+KTA1+KTB1+KTA2+KTB2+KTA3+KTB3

	doeclimvars.Adoe[1,1] = 1 - deltat/(2.*doeclimvars.taucfl) - deltat/(2.*doeclimvars.taukls)
    doeclimvars.Adoe[1,2] =  deltat/(2.*doeclimvars.taukls)*bsi
    doeclimvars.Adoe[2,1] =  deltat/(2.*doeclimvars.tauksl)
    doeclimvars.Adoe[2,2] = 1 - deltat/(2.*doeclimvars.taucfs) - deltat/(2.*doeclimvars.tauksl)*bsi + doeclimvars.Ker[nsteps]*fso*sqrt(deltat/doeclimvars.taudif)    
    doeclimvars.Adoe=doeclimvars.Adoe+Cdoe
end

#!------------------------------------------------------------------------------
#
#!------------------------------------------------------------------------------
#SUBROUTINE doeclimtimestep_simple(n,forcing,temp)
#!  ==========================================================================
#! | Simple climate model DOECLIM
#! |
#! | calculates sea surface and land air temperature response to radiative 
#! | forcing based on an energy balance model with 1-D diffusion ocean
#! |
#! | *** computes single time step ***
#! | *** initialize with init_doeclim ***
#! | *** then iterate this function ***
#! |
#! | Input:
#! |       n:        current time step
#! |       forcing:  global radiative forcing (top of atmosphere) (W/m^2)
#! |
#! | Output:
#! |       temp: global mean temperature anomaly (K), relative to preindustrial
#! |
#! | Assumptions: 
#! |       land surface temperature = land air temperature
#! |       mixed layer temperature  = sea surface temperatures 
#! |                                = marine air temperature divided by bsi
#!  ==========================================================================#
#
#    implicit none#
#
#    integer(i4b),intent(IN) :: n
#    real(DP),intent(OUT) :: temp
#    real(DP), dimension(2) :: DQ, DPAST, QC, DTEAUX
#    real(DP), dimension(nsteps) :: QL, Q0
#    real(DP), dimension(2,nsteps) :: DTE
#    real(DP) :: DelQL, DelQ0
#    real(DP) :: forcing
#    integer(i4b) :: i
#
#     DTE(1,:) = temp_landair
#     DTE(2,:) = temp_sst
#!% assume land and ocean forcings are equal to global forcing
#     QL = forcing !-forcing(1)!-1.4D-1!-FORC(1)    !forcing
#     Q0 = forcing !-forcing(1)!-1.4D-1!-FORC(1)    !forcing
#
#   if (n.gt.2) then
#     DelQL = QL(n) - QL(n-1)
#     DelQ0 = Q0(n) - Q0(n-1)
#
#!    % Assumption: linear forcing change between n and n+1
#!     do i=1,nsteps
#     QC(1) = (DelQL/cal*(1./taucfl+1./taukls)-bsi*DelQ0/cas/taukls)
#     QC(2) = (DelQ0/cas*(1./taucfs+bsi/tauksl)-DelQL/cal/tauksl)
#!     enddo
#     QC = QC* deltat^2/12.
#!    % -------------------------- INITIAL CONDITIONS ------------------------
#!    % Initialization of temperature and forcing vector:
#!    %Factor 1/2 in front of Q in Equation A.27, EK05, and Equation 2.3.27, TK07 is a typo!
#!    %Assumption: linear forcing change between n and n+1
#     DQ(1) = 0.5d0*deltat/cal*(QL(n)+QL(n-1))
#     DQ(2) = 0.5d0*deltat/cas*(Q0(n)+Q0(n-1))
#     DQ = DQ + QC
#
#!    % -------------- SOLVE MODEL ------------------------------------
#!    % Calculate temperatures
#!     DPAST = zeros(2,1)
#     DPAST = 0.0d0
#     do i=1,n-1
#        DPAST(2) = DPAST(2)+DTE(2,i)*Ker(nsteps-n+i)
#     enddo
#     DPAST(2) = DPAST(2)*fso * sqrt(deltat/taudif)
#!     DPAST(2) = fso * sqrt(deltat/taudif) * DTE(2,(1:n-1))*Ker(nsteps-n+1:nsteps-1)  !was transposed
#
#!     DTE(:,n) = IB * ( DQ + DPAST + A*DTE(:,n-1) )
#
#     DTEAUX(1) = Adoe(1,1)*DTE(1,n-1)+Adoe(1,2)*DTE(2,n-1)
#     DTEAUX(2) = Adoe(2,1)*DTE(1,n-1)+Adoe(2,2)*DTE(2,n-1)
#
#     DTE(1,n) = IB(1,1)*(DQ(1)+DPAST(1)+DTEAUX(1))+                  &
#                IB(1,2)*(DQ(2)+DPAST(2)+DTEAUX(2))
#     DTE(2,n) = IB(2,1)*(DQ(1)+DPAST(1)+DTEAUX(1))+                  &
#                IB(2,2)*(DQ(2)+DPAST(2)+DTEAUX(2))
#
#   else
#!% handle initial values
#
#    DTE(:,1) = 0D0
#!    DTE(:,1) = T_ATM_0!0D0
#!     DTE(:,1) = 0.6024D0
#
#     DelQL = QL(2) - QL(1)
#     DelQ0 = Q0(2) - Q0(1)
#
#     QC(1) = DelQL/cal*(1./taucfl+1./taukls)-bsi*DelQ0/cas/taukls
#     QC(2) = DelQ0/cas*(1./taucfs+bsi/tauksl)-DelQL/cal/tauksl
#     QC = QC*deltat**2/12.
#
#     DQ(1) = 0.5d0*deltat/cal*(QL(2)+QL(1))
#     DQ(2) = 0.5d0*deltat/cas*(Q0(2)+Q0(1))
#     DQ= DQ + QC
#
#!     DTE(:,2) = IB * ( DQ +  A*DTE(:,1) )
#     DTEAUX(1)=Adoe(1,1)*DTE(1,1)+Adoe(1,2)*DTE(2,1)
#     DTEAUX(2)=Adoe(2,1)*DTE(1,1)+Adoe(2,2)*DTE(2,1)


#     DTE(1,2) = IB(1,1)*(DQ(1) +DTEAUX(1))+IB(1,2)*(DQ(2)+DTEAUX(2))
#     DTE(2,2) = IB(2,1)*(DQ(1) +DTEAUX(1))+IB(2,2)*(DQ(2)+DTEAUX(2))
#   endif

#     temp_landair(n) = DTE(1,n)
#     temp_sst(n) = DTE(2,n)
#     temp = flnd*temp_landair(n) + (1.-flnd)*bsi*temp_sst(n)

#!% Calculate ocean heat uptake [W/m^2]
#!% heatflux(n) captures in the heat flux in the period between n-1 and n
#!% Numerical implementation of Equation 2.7, EK05, or Equation 2.3.13, TK07)
#!% ------------------------------------------------------------------------


#   if (n.gt.1) then
#     heatflux_mixed(n) = cas*( DTE(2,n)-DTE(2,n-1) )

#!     heatflux_interior(n) = cas*fso/sqrt(taudif*deltat)*            &
#!     ( 2.*DTE(2,n) - DTE(2,(1:n-1))*Ker(nsteps-n+2:nsteps) )  !was transposed
#!     heatflux_interior(1)=0.
#     do i=1,n-1
#        heatflux_interior(n) = heatflux_interior(n)+DTE(2,i)*Ker(nsteps-n+1+i)
#     enddo
#     heatflux_interior(n) = cas*fso/sqrt(taudif*deltat)*(2.*DTE(2,n) -       &
#                            heatflux_interior(n))

#     heat_mixed(n) = heat_mixed(n-1) +heatflux_mixed(n) *(powtoheat*deltat)

#     heat_interior(n) = heat_interior(n-1) + heatflux_interior(n) *      &
#                        (fso*powtoheat*deltat)

#   else
#!% handle initial values
#     heatflux_mixed(1) = 0.0d0
#     heatflux_interior(1) = 0.0d0

#    heat_mixed(1) = 0.0d0
#     heat_interior(1) = 0.0d0
#   endif

#END SUBROUTINE doeclimtimestep_simple
end