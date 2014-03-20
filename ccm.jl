#!  CCM:   Carbon Cycle Model
#!
#!  Copyright (C) 2009 D. Ricciuto, B. Tuttle, K. Keller
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
#! Nonlinear impulse response carbon/climate model
#! based on the NICCS model from Hooss et al. (2001)
#! original version by DMRicciuto 7/16/2004
#!
#! Notes from model-code-draft-3/src/model.f90 (D. McInerney):
#!   Terrestrial model:  4 box model based on Meyer et al. (1999)
#!   Model to be calibrated based on flux tower data synthesis
#!   Oceanic model currently using linear IRFs
#!   Climate model calibrated to Hamburg AOGCM (T only)
#!
#! See model details in:
#!   Ricciuto, D. M., K. J. Davis, and K. Keller (2008), A Bayesian calibration 
#!       of a simple carbon cycle model: The role of observations in estimating 
#!       and reducing uncertainty, Global Biogeochem. Cycles, 22, GB2030, 
#!       doi:10.1029/2006GB002908.
#!
#!------------------------------------------------------------------------------
#!  13 Mar 2009  Brian Tuttle <btuttle@psu.edu> received 
#!               globalinversion_fortran/model.f from Dan Ricciuto.
#!   May 2009    Rewrote model() subroutine as CCM.f90 module, including 
#!               initialization, (de)allocation, and CC_model subroutines.
#!   Aug 2009    Incorporated CCM.f90 into EarthSystem module.
#!               Removed usetemp logical switch as well as the simple impulse
#!               response carbon/climate model in lieue of externally 
#!               computed temperature forcing.
#!------------------------------------------------------------------------------

module ccm

# Define factors used in the calculation.
const r3f = 45.0 / 120.0
const tp1f = 35.0 / 60.0
const tp2f = 25.0 / 60.0
const tp3f = 1.0 / 12.0

# Define model parameters.
const hs = 64.0    #layer depths (m)
const h1 = 672.0
const h2 = 419.0
const h3 = 1136.0
const h4 = 2382.0
const n3 = 9.04 
const n4 = 6.32
const npp0 = 60.0             # [GtC/yr]

type ccmpar
    nsteps::Int
    deltat::Float64
    Clim_sens::Float64
    Q10::Float64
    Beta::Float64
    Eta::Float64 #default 16.88,   diffusion coeffs [m/yr]
    temp::Vector{Float64}
    CO2_emissions::Vector{Float64}
    anomtable::Array{Float64,2}
end

type ccmvar
    tpools::Array{Float64,2}
    ocanom::Array{Float64,2}
    atmco2::Vector{Float64}
    landflux::Vector{Float64}
    atm_oc_flux::Vector{Float64}

    function ccmvar(p::ccmpar)
        vars = new(
            zeros(p.nsteps+1,4),
            zeros(p.nsteps+1,4),
            zeros(p.nsteps+1),
            zeros(p.nsteps),
            zeros(p.nsteps))
        return vars
    end
end

#! Input parameters:
#    real(DP) :: Cs          ! Climate sensitivity
#    real(DP) :: Q10         ! Respiration Temperature sens.
#    real(DP) :: Beta        ! Carbon Fertilization param.
#    real(DP) :: Eta         ! Thermocline transfer velocity

#    real(DP), dimension(2) :: trm
#    real(DP), dimension(2) :: err
#    real(DP), dimension(:,:), allocatable, public :: anomtable
#    integer(i4b), dimension(2), public :: ATsize 

#! Model parameters:
#    real(DP) :: n2, n3, n4, hs, h1, h2, h3, h4, npp0
#    real(DP), dimension(:,:), allocatable :: tpools
#    real(DP), dimension(:,:), allocatable :: ocanom
#    real(DP), dimension(4) :: Ftp, Goc

#! Model factors:
#!    real(DP) :: deltat

#! Climate model variables:
#!    real(DP) :: a1, a2, tao1, tao2
#    real(DP), dimension(:), allocatable, public :: atmco2
#    real(DP), dimension(:), allocatable, public :: landflux
#    real(DP), dimension(:), allocatable, public :: atm_oc_flux
#    real(DP), dimension(:), allocatable, public :: emissions    ![GtC/yr]

#    public :: init_CCM_arrays, init_CCM_parameters, CC_model, dealloc_CCM
#    public :: alloc_anomtab, dealloc_anomtab

function init(p::ccmpar, v::ccmvar)
    #  =========================================================================
    # |  Allocate and initialize Carbon Cycle Model arrays.  Set first array    |
    # |  elements to preindustrial values.                                      |
    #  =========================================================================


    #! Initialize terrestrial pools (equilibrium preindustrial values)
    #    tpools = 0.0
    #    Ftp = 0.0

    v.tpools[1,1] = 100.0     # Non-woody vegetation [GtC]
    v.tpools[1,2] = 500.0     # Woody vegetation [GtC]
    v.tpools[1,3] = 120.0     # Detritus [GtC]
    v.tpools[1,4] = 1500.0    # Soil carbon [GtC]

    #    landflux = 0.0
      
    #    ocanom = 0.0
    #    Goc = 0.0

    #    atmco2(:) = 0.0
    v.atmco2[1] = 285.2        # [ppm]

    #! Assign inputs to global variables.
    #Cs = Clim_sens
    #Q10 = Soil_resp
    #Beta = Carb_fert
    #Eta = Therm_diff
end

function timestep(p::ccmpar, v::ccmvar, t::Int)
#  =========================================================================
# |  Carbon Cycle Model, single time step
# |
# |  Input parameters:
# |     t:      time index
# |     temp:   temperature forcing [K]
# |     CO2_emissions:  [GtC/yr]
# |
# |  Optional output:
# |     diagnostic: array of model variables
#  =========================================================================
      
#      real(DP) :: fracinoc
     
    Q10temp = p.Q10^(p.temp[t]/10.0)

    # Calculate Net Primary Productivity.   (eq2, Ricciuto 2008)
    npp = npp0 * (1.0+p.Beta*log(v.atmco2[t]/v.atmco2[1]))

    # Calculate Heterotrophic respiration rate.     (eq3, Ricciuto 2008)
    resp3 = v.tpools[t,3] * r3f * Q10temp
    resp4 = v.tpools[t,4] * 0.01 * Q10temp
    resp_h = resp3+resp4

    v.landflux[t] = resp_h - npp
	
    # Set terrestrial pool sizes for next timestep	
    Ftp = Array(Float64,4)
    Ftp[1] = npp*tp1f - 0.35*v.tpools[t,1]
    Ftp[2] = npp*tp2f - 0.05*v.tpools[t,2]
    Ftp[3] = 0.35*v.tpools[t,1] + 0.04*v.tpools[t,2] - tp3f*v.tpools[t,3] - resp3
    Ftp[4] = 0.01*v.tpools[t,2] + tp3f*v.tpools[t,3] - resp4
      
    v.tpools[t+1,:] = vec(v.tpools[t,:]) + p.deltat*Ftp[:]

    netemissions = p.CO2_emissions[t] + v.landflux[t]
	   

    # Find fracinoc using the ocean anomaly table.
    fracinoc = anom_interp(p.anomtable, p.temp[t], v.ocanom[t,1]+netemissions)
    
    # Compute carbon anomaly in ocean mixed layer/atmosphere
    Goc = Array(Float64,4)
    Goc[1] = netemissions-(p.Eta/hs)*v.ocanom[t,1]* fracinoc+(p.Eta/h2)*v.ocanom[t,2]
    # Compute carbon anomaly in ocean layers 2-4
    Goc[2] = (p.Eta/hs)*v.ocanom[t,1]*fracinoc- ((p.Eta+n3)/h2)*v.ocanom[t,2]+(n3/h3)*v.ocanom[t,3]
    Goc[3] = (n3/h2)*v.ocanom[t,2]-((n3+n4)/h3)*v.ocanom[t,3]+ (n4/h4)*v.ocanom[t,4]
    Goc[4] = (n4/h3)*v.ocanom[t,3] - (n4/h4)* v.ocanom[t,4]

    v.ocanom[t+1,:] = vec(v.ocanom[t,:]) + p.deltat*Goc[:]
    
    # Compute flux into ocean
    v.atm_oc_flux[t] = -((v.ocanom[t+1,1]-v.ocanom[t,1])*fracinoc + v.ocanom[t+1,2]-v.ocanom[t,2]+v.ocanom[t+1,3]-v.ocanom[t,3] + v.ocanom[t+1,4]-v.ocanom[t,4]) / p.deltat
    v.atmco2[t+1] = v.atmco2[1]+(v.ocanom[t+1,1]*(1-fracinoc))/2.13

    #emissions(t) = CO2_emissions
end

function anom_interp(anomtable, ref_temp::Float64, ref_emis::Float64)
    #  =========================================================================
    # |  A simple linear 2D interpolation over the ocean anomaly table to find  |
    # |  the fracinoc variable.                                                 |
    #  =========================================================================

    #integer(i4b) :: templox, temphix, emislox, emishix 
    #real(DP) :: tempx, emisx, FIO_tlo, FIO_thi
    ATsize = size(anomtable)

    # Upper and lower bound indices (integers).
     templox::Int = min(max(floor(ref_temp*10.0+10.0)+1,1),ATsize[1])
     temphix::Int = max(min(iceil(ref_temp*10.0+10.0)+1,ATsize[1]),1)
     emislox::Int = min(max(floor((ref_emis)/2.0)+1,1),ATsize[2])
     emishix::Int = max(min(iceil((ref_emis)/2.0)+1,ATsize[2]),1)
    # Target indices (reals).
     tempx::Float64 = (ref_temp*10.0+10.0)+1.0
     emisx::Float64 = ((ref_emis)/2.0)+1.0

    # First interpolate anomtable in the emission direction.
    if emislox == emishix
        FIO_tlo = anomtable[templox,emislox]
        FIO_thi = anomtable[temphix,emislox]
    else
        FIO_tlo = (anomtable[templox,emishix]-anomtable[templox,emislox]) * (emisx-emislox) / real(emishix-emislox) + anomtable[templox,emislox]
        FIO_thi = (anomtable[temphix,emishix]-anomtable[temphix,emislox]) * (emisx-emislox) / real(emishix-emislox) + anomtable[temphix,emislox]
    end

    # Then interpolate the result in the temperature direction.
    if templox == temphix then
        frac_in_ocean = FIO_tlo
    else
        frac_in_ocean = (FIO_thi - FIO_tlo) * (tempx-templox) / real(temphix-templox) + FIO_tlo
    end

    return frac_in_ocean
end

end

xpar = ccm.ccmpar(200,1,3,1.126598,0.2916273,172.2809,zeros(200),zeros(200),zeros(100,10000))
xpar.temp[:] = .7
xpar.CO2_emissions[:] = 7
xpar.anomtable[:,:] = 0.0


xvar = ccm.ccmvar(xpar)
ccm.init(xpar, xvar)

for t=1:200    
    ccm.timestep(xpar, xvar, t)
end