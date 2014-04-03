#  CCM:   Carbon Cycle Model
#------------------------------------------------------------------------------
# Nonlinear impulse response carbon/climate model
# based on the NICCS model from Hooss et al. (2001)
# original version by DMRicciuto 7/16/2004
#
# Notes from model-code-draft-3/src/model.f90 (D. McInerney):
#   Terrestrial model:  4 box model based on Meyer et al. (1999)
#   Model to be calibrated based on flux tower data synthesis
#   Oceanic model currently using linear IRFs
#   Climate model calibrated to Hamburg AOGCM (T only)
#
# See model details in:
#   Ricciuto, D. M., K. J. Davis, and K. Keller (2008), A Bayesian calibration 
#       of a simple carbon cycle model: The role of observations in estimating 
#       and reducing uncertainty, Global Biogeochem. Cycles, 22, GB2030, 
#       doi:10.1029/2006GB002908.
#
#------------------------------------------------------------------------------
#  13 Mar 2009  Brian Tuttle <btuttle@psu.edu> received 
#               globalinversion_fortran/model.f from Dan Ricciuto.
#   May 2009    Rewrote model() subroutine as CCM.f90 module, including 
#               initialization, (de)allocation, and CC_model subroutines.
#   Aug 2009    Incorporated CCM.f90 into EarthSystem module.
#               Removed usetemp logical switch as well as the simple impulse
#               response carbon/climate model in lieue of externally 
#               computed temperature forcing.
#------------------------------------------------------------------------------

module ccmcomponent
using IAMF

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

@defcomp ccm begin
    addParameter(deltat,Float64)
    addParameter(Clim_sens,Float64)
    addParameter(Q10,Float64)
    addParameter(Beta,Float64)
    addParameter(Eta,Float64) #default 16.88,   diffusion coeffs [m/yr]
    addParameter(temp,Float64, index=[time])
    addParameter(CO2_emissions,Float64, index=[time])
    addParameter(anomtable,Array{Float64,2})

    addVariable(tpools, Float64, index=[time,4])
    addVariable(ocanom, Float64, index=[time,4])
    addVariable(atmco2, Float64, index=[time])
    addVariable(landflux, Float64, index=[time])
    addVariable(atm_oc_flux, Float64, index=[time])
    addVariable(Ftp, Float64, index=[4])
    addVariable(Goc, Float64, index=[4])
end

function init(s::ccm)
    p = s.Parameters
    v = s.Variables
    v.tpools[1,1] = 100.0     # Non-woody vegetation [GtC]
    v.tpools[1,2] = 500.0     # Woody vegetation [GtC]
    v.tpools[1,3] = 120.0     # Detritus [GtC]
    v.tpools[1,4] = 1500.0    # Soil carbon [GtC]

    v.atmco2[1] = 285.2        # [ppm]
end

function timestep(s::ccm, t::Int)
    p = s.Parameters
    v = s.Variables

    Q10temp = p.Q10^(p.temp[t]/10.0)

    # Calculate Net Primary Productivity.   (eq2, Ricciuto 2008)
    npp = npp0 * (1.0+p.Beta*log(v.atmco2[t]/v.atmco2[1]))

    # Calculate Heterotrophic respiration rate.     (eq3, Ricciuto 2008)
    resp3 = v.tpools[t,3] * r3f * Q10temp
    resp4 = v.tpools[t,4] * 0.01 * Q10temp
    resp_h = resp3+resp4

    v.landflux[t] = resp_h - npp
    
    # Set terrestrial pool sizes for next timestep  
    v.Ftp[1] = npp*tp1f - 0.35*v.tpools[t,1]
    v.Ftp[2] = npp*tp2f - 0.05*v.tpools[t,2]
    v.Ftp[3] = 0.35*v.tpools[t,1] + 0.04*v.tpools[t,2] - tp3f*v.tpools[t,3] - resp3
    v.Ftp[4] = 0.01*v.tpools[t,2] + tp3f*v.tpools[t,3] - resp4
      
    if t<s.nsteps
        for i=1:4
            v.tpools[t+1,i] = v.tpools[t,i] + p.deltat*v.Ftp[i]
        end
    end

    netemissions = p.CO2_emissions[t] + v.landflux[t]      

    # Find fracinoc using the ocean anomaly table.
    fracinoc = anom_interp(p.anomtable, p.temp[t], v.ocanom[t,1]+netemissions)
    
    # Compute carbon anomaly in ocean mixed layer/atmosphere
    v.Goc[1] = netemissions-(p.Eta/hs)*v.ocanom[t,1]* fracinoc+(p.Eta/h2)*v.ocanom[t,2]
    # Compute carbon anomaly in ocean layers 2-4
    v.Goc[2] = (p.Eta/hs)*v.ocanom[t,1]*fracinoc- ((p.Eta+n3)/h2)*v.ocanom[t,2]+(n3/h3)*v.ocanom[t,3]
    v.Goc[3] = (n3/h2)*v.ocanom[t,2]-((n3+n4)/h3)*v.ocanom[t,3]+ (n4/h4)*v.ocanom[t,4]
    v.Goc[4] = (n4/h3)*v.ocanom[t,3] - (n4/h4)* v.ocanom[t,4]

    if t<s.nsteps
        for i=1:4
            v.ocanom[t+1,i] = v.ocanom[t,i] + p.deltat*v.Goc[i]
        end
    
        # Compute flux into ocean
        v.atm_oc_flux[t] = -((v.ocanom[t+1,1]-v.ocanom[t,1])*fracinoc + v.ocanom[t+1,2]-v.ocanom[t,2]+v.ocanom[t+1,3]-v.ocanom[t,3] + v.ocanom[t+1,4]-v.ocanom[t,4]) / p.deltat

        v.atmco2[t+1] = v.atmco2[1]+(v.ocanom[t+1,1]*(1-fracinoc))/2.13
    end
end

function anom_interp(anomtable, ref_temp::Float64, ref_emis::Float64)
    # A simple linear 2D interpolation over the ocean anomaly table to find
    # the fracinoc variable.
    
    ATsize = size(anomtable)

    # Upper and lower bound indices (integers).
    templox::Int = min(max(ifloor(ref_temp*10.0+10.0)+1,1),ATsize[1])
    temphix::Int = max(min(iceil(ref_temp*10.0+10.0)+1,ATsize[1]),1)
    emislox::Int = min(max(ifloor((ref_emis)/2.0)+1,1),ATsize[2])
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