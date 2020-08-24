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
# of climate change. Ph.D. thesiv. University of Potsdam, 256 pp.
# opuv.kobv.de/ubp/volltexte/2005/561/
#
# Model equations are described in EK05 and (Reference TK07):
# Tanaka, K, Kriegler, E, Bruckner, T, Hooss, G, Knorr, W, Raddatz, T (2007)
# Aggregated carbon cycle, atmospheric chemistry, and climate model (ACC2):
# Description of the forward and inverse modes, Reports on Earth System Science ! 40/2007,
# Max Planck Institute for Meteorology, Hamburg, 199 pp.
# www.mpimet.mpg.de/fileadmin/publikationen/Reports/BzE_40.pdf
#
# ==============================================================================
#
# Updates:
# 22.05.2007 Hammer-Hollingsworth numerical correction included (EK)
# 23.05.2007 Ocean heat uptake added (EK)
# 12.02.2008 Translated to Fortran90 (Marlos Goes <mpg14@psu.edu>)
# 15.08.2009 Written as Fortran90 module (Brian Tuttle <btuttle@psu.edu>)
#
# ==============================================================================
#
# Global Parameters:
#   ak      slope coeff. for land-sea heat exchange
#   bk      interv. coeff. for land-sea heat exch.
#   bsi     marine air warming enhancement
#   cal     heat cap. of land-troposph. system
#   cas     heat cap. of ocean ML-troposph.
#   csw     specific heat capacity of 1m^3 seawater [Wa/m^3/K]
#   deltat  time step size [years]
#   flnd    land fraction
#   fso     ocean frac. area below 60m
#   kcon    conversion factor [cm2/s->m2/a]
#   F2x_CO₂ forcing increase from CO2 doubling [W/m^2]
#   rlam    clim senv. over land enhancement
#   zbot    depth of interior ocean
#
#   temp_landair:       land air temperature anomaly (K)
#   temp_sst:           sea surface temperature anomaly (K)
#   heat_mixed:         mixed layer heat anomaly (10^22 J)
#   heat_interior:      interior ocean heat anomaly (10^22 J)
#   heatflux_mixed:     heat uptake of the mixed layer (W/m^2)
#   heatflux_interior:  heat uptake of the interior ocean (W/m^2)
#
# ==============================================================================

const ak   = 0.31
const bk   = 1.59
const bsi  = 1.3
const cal  = 0.52
const cas  = 7.80
const csw  = 0.13
const flnd = 0.29
const fso  = 0.95
const kcon = 3155.8
const rlam = 1.43
const zbot = 4000
const earth_area = 5100656e8      # [m^2]
const secs_per_Year = 31556926.0

@defcomp doeclim begin
    deltat = Parameter()
    # climate sensitivity to 2xCO2 (K); default = 3
    t2co = Parameter()
    # vertical ocean diffusivity (cm^2 s^-1); default = 0.55
    kappa = Parameter()
    forcing = Parameter(index=[time])
    F2x_CO₂ = Parameter(default=3.7)

    temp = Variable(index=[time])
    temp_landair = Variable(index=[time])
    temp_sst = Variable(index=[time])
    heat_mixed = Variable(index=[time])
    heat_interior = Variable(index=[time])
    heatflux_mixed = Variable(index=[time])
    heatflux_interior = Variable(index=[time])
    Ker = Variable(index=[time])

    taucfl = Variable()
    taukls = Variable()
    taucfs = Variable()
    tauksl = Variable()
    taudif = Variable()
    taubot = Variable()
    powtoheat = Variable()

    IB = Variable(index=[2,2])
    Adoe = Variable(index=[2,2])

    KT0 = Variable(index=[time])
    KTA1 = Variable(index=[time])
    KTB1 = Variable(index=[time])
    KTA2 = Variable(index=[time])
    KTB2 = Variable(index=[time])
    KTA3 = Variable(index=[time])
    KTB3 = Variable(index=[time])

    Cdoe = Variable(index=[2,2])
    Baux = Variable(index=[2,2])

    function init(p, v, d)   
        # DEPENDENT MODEL PARAMETERS
        ocean_area = (1.0-flnd)*earth_area
    
        v.powtoheat = ocean_area*secs_per_Year / 1e22
    
        cnum = rlam*flnd + bsi * (1.0-flnd)
    
        cden = rlam * flnd - ak *(rlam-bsi)
    
        # vertical diffusivity in [m^2/a]
    
        keff = kcon * p.kappa
    
        # climate feedback strength over land
    
        cfl = flnd *cnum/cden*p.F2x_CO₂/p.t2co-bk*(rlam-bsi)/cden
    
        # climate feedback strength over ocean
    
        cfs = (rlam * flnd - ak / (1-flnd) * (rlam-bsi)) * cnum / cden * p.F2x_CO₂ / p.t2co + rlam * flnd / (1-flnd) * bk * (rlam - bsi) / cden
    
        # land-sea heat exchange coefficient
    
        kls = bk * rlam * flnd / cden - ak * flnd * cnum / cden * p.F2x_CO₂ / p.t2co
    
        # interior ocean warming time scale
    
        v.taubot = zbot^2 / keff
    
        # ocean heat diff. time scale
    
        v.taudif = cas^2 / csw^2 * pi / keff
    
        # ocean response time scale
    
        v.taucfs = cas / cfs
    
        # land response time scale
    
        v.taucfl = cal / cfl
    
        # sea-land heat exchange time scale
    
        v.tauksl  = (1-flnd) * cas / kls
    
        # land-sea heat exchange time scale
    
        v.taukls  = flnd * cal / kls
    
        # Zeroth Order
    
        v.KT0[end] = 4-2*sqrt(2.)
    
        # First Order
    
        v.KTA1[end] = -8*exp(-v.taubot/p.deltat) + 4*sqrt(2.)*exp(-0.5*v.taubot/p.deltat)
    
        v.KTB1[end] = 4*sqrt(pi*v.taubot/p.deltat) * (1+erf(sqrt(0.5*v.taubot/p.deltat)) - 2*erf(sqrt(v.taubot/p.deltat)))
    
        # Second order
    
        v.KTA2[end] =  8*exp(-4*v.taubot/p.deltat) - 4*sqrt(2.)*exp(-2*v.taubot/p.deltat)
    
        v.KTB2[end] = -8*sqrt(pi*v.taubot/p.deltat) * (1+erf(sqrt(2*v.taubot/p.deltat)) - 2*erf(2*sqrt(v.taubot/p.deltat)) )
    
        # Third Order
    
        v.KTA3[end] = -8*exp(-9*v.taubot/p.deltat) + 4*sqrt(2.)*exp(-4.5*v.taubot/p.deltat)
    
        v.KTB3[end] = 12*sqrt(pi*v.taubot/p.deltat) * (1 +erf(sqrt(4.5*v.taubot/p.deltat)) - 2*erf(3*sqrt(v.taubot/p.deltat)) )
    
        # Hammer and Hollingsworth correction (Equation 2.3.27, TK07):
        # Switched on (To switch off, comment out lines below)
        v.Cdoe[1,1] = (1/v.taucfl^2+1/v.taukls^2+2/v.taucfl/v.taukls+bsi/v.taukls/v.tauksl) * (p.deltat^2/12.)
        v.Cdoe[1,2] = (-bsi/v.taukls^2-bsi/v.taucfl/v.taukls-bsi/v.taucfs/v.taukls-bsi^2/v.taukls/v.tauksl) * (p.deltat^2/12.)
        v.Cdoe[2,1] = (-bsi/v.tauksl^2-1/v.taucfs/v.tauksl-1/v.taucfl/v.tauksl-1/v.taukls/v.tauksl) * (p.deltat^2/12.)
        v.Cdoe[2,2] = (1/v.taucfs^2+bsi^2/v.tauksl^2+2*bsi/v.taucfs/v.tauksl+bsi/v.taukls/v.tauksl) * (p.deltat^2/12.)
    
        # Matrices of difference equation system B*T(i+1) = Q(i) + A*T(i)
        # T = (TL,TO)
        # (Equation A.27, EK05, or Equations 2.3.24 and 2.3.27, TK07)
        v.Baux[1,1] = 1. + p.deltat/(2*v.taucfl) + p.deltat/(2*v.taukls)+v.Cdoe[1,1]
        v.Baux[1,2] = -p.deltat/(2*v.taukls)*bsi+v.Cdoe[1,2]
        v.Baux[2,1] = -p.deltat/(2*v.tauksl)+v.Cdoe[2,1]
        v.Baux[2,2] = 1. + p.deltat/(2*v.taucfs) + p.deltat/(2*v.tauksl)*bsi + 2*fso*sqrt(p.deltat/v.taudif)+v.Cdoe[2,2]
    
        v.IB[:] = inv(v.Baux)

        # TODO find better solution
        nsteps = length(v.Ker)
    
        for i = 1:nsteps-1
    
            # Zeroth Order
    
            v.KT0[TimestepIndex(i)] = 4*sqrt(Float64(nsteps+1-i)) - 2*sqrt(Float64(nsteps+2-i)) - 2*sqrt(Float64(nsteps-i))
    
            # First Order
    
            v.KTA1[TimestepIndex(i)] = -8*sqrt(Float64(nsteps+1-i)) * exp(-v.taubot/p.deltat/(nsteps+1-i)) + 4*sqrt(Float64(nsteps+2-i)) *exp(-v.taubot/p.deltat/(nsteps+2-i)) + 4*sqrt(Float64(nsteps-i)) *exp(-v.taubot/p.deltat/(nsteps-i))
    
            v.KTB1[TimestepIndex(i)] = 4*sqrt(pi*v.taubot/p.deltat) * (erf(sqrt(v.taubot/p.deltat/(nsteps-i))) + erf(sqrt(v.taubot/p.deltat/(nsteps+2-i))) - 2*erf(sqrt(v.taubot/p.deltat/(nsteps+1-i))) )
    
            # Second Order
    
            v.KTA2[TimestepIndex(i)] =  8*sqrt(Float64(nsteps+1-i)) * exp(-4*v.taubot/p.deltat/(nsteps+1-i))- 4*sqrt(Float64(nsteps+2-i))*exp(-4*v.taubot/p.deltat/(nsteps+2-i))- 4*sqrt(Float64(nsteps-i)) * exp(-4*v.taubot/p.deltat/(nsteps-i))
    
            v.KTB2[TimestepIndex(i)] = -8*sqrt(pi*v.taubot/p.deltat) * (erf(2*sqrt(v.taubot/p.deltat/(Float64(nsteps-i)))) + erf(2*sqrt(v.taubot/p.deltat/Float64(nsteps+2-i))) -       2*erf(2*sqrt(v.taubot/p.deltat/Float64(nsteps+1-i))) )
    
            # Third Order
    
            v.KTA3[TimestepIndex(i)] = -8*sqrt(Float64(nsteps+1-i)) *exp(-9*v.taubot/p.deltat/(nsteps+1-i)) + 4*sqrt(Float64(nsteps+2-i))*exp(-9*v.taubot/p.deltat/(nsteps+2-i)) + 4*sqrt(Float64(nsteps-i))*exp(-9*v.taubot/p.deltat/(nsteps-i))
    
            v.KTB3[TimestepIndex(i)] = 12*sqrt(pi*v.taubot/p.deltat) * (erf(3*sqrt(v.taubot/p.deltat/(nsteps-i))) + erf(3*sqrt(v.taubot/p.deltat/(nsteps+2-i))) - 2*erf(3*sqrt(v.taubot/p.deltat/(nsteps+1-i))) )
        end
    
        for i=1:nsteps
            v.Ker[TimestepIndex(i)] = v.KT0[TimestepIndex(i)]+v.KTA1[TimestepIndex(i)]+v.KTB1[TimestepIndex(i)]+v.KTA2[TimestepIndex(i)]+v.KTB2[TimestepIndex(i)]+v.KTA3[TimestepIndex(i)]+v.KTB3[TimestepIndex(i)]
        end

        v.Adoe[1,1] = 1 - p.deltat/(2*v.taucfl) - p.deltat/(2*v.taukls) + v.Cdoe[1,1]
        v.Adoe[1,2] =  p.deltat/(2*v.taukls)*bsi + v.Cdoe[1,2]
        v.Adoe[2,1] =  p.deltat/(2*v.tauksl) + v.Cdoe[2,1]
        v.Adoe[2,2] = 1 - p.deltat/(2*v.taucfs) - p.deltat/(2*v.tauksl)*bsi + v.Ker[TimestepIndex(nsteps)]*fso*sqrt(p.deltat/v.taudif) + v.Cdoe[2,2]
    end
    
    function run_timestep(p, v, d, n)
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
        DTE1 = v.temp_landair
        DTE2 = v.temp_sst

        # TODO find better solution
        nsteps = length(v.Ker)
    
        # assume land and ocean forcings are equal to global forcing
        QL = p.forcing
        Q0 = p.forcing

        if !is_first(n)
            DelQL = QL[n] - QL[n-1]
            DelQ0 = Q0[n] - Q0[n-1]
    
            # Assumption: linear forcing change between n and n+1
            QC1 = (DelQL/cal*(1/v.taucfl+1/v.taukls)-bsi*DelQ0/cas/v.taukls)
            QC2 = (DelQ0/cas*(1/v.taucfs+bsi/v.tauksl)-DelQL/cal/v.tauksl)
    
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

            # TODO convert this to proper timestep handling
            for i=1:n.t-1
                DPAST2 = DPAST2+DTE2[TimestepIndex(i)]*v.Ker[TimestepIndex(nsteps-n.t+i)]
            end
    
            DPAST2 = DPAST2*fso * sqrt(p.deltat/v.taudif)
    
            DTEAUX1 = v.Adoe[1,1]*DTE1[n-1]+v.Adoe[1,2]*DTE2[n-1]
            DTEAUX2 = v.Adoe[2,1]*DTE1[n-1]+v.Adoe[2,2]*DTE2[n-1]
    
            DTE1[n] = v.IB[1,1]*(DQ1+DPAST1+DTEAUX1)+v.IB[1,2]*(DQ2+DPAST2+DTEAUX2)
            DTE2[n] = v.IB[2,1]*(DQ1+DPAST1+DTEAUX1)+v.IB[2,2]*(DQ2+DPAST2+DTEAUX2)
    
    
            # Calculate ocean heat uptake [W/m^2]
            # heatflux(n) captures in the heat flux in the period between n-1 and n
            # Numerical implementation of Equation 2.7, EK05, or Equation 2.3.13, TK07)
            # ------------------------------------------------------------------------
            v.heatflux_mixed[n] = cas*( DTE2[n]-DTE2[n-1] )
    
            v.heatflux_interior[n] = 0.
            for i=1:n.t-1
                v.heatflux_interior[n] = v.heatflux_interior[n]+DTE2[TimestepIndex(i)]*v.Ker[TimestepIndex(nsteps-n.t+1+i)]
            end
            v.heatflux_interior[n] = cas*fso/sqrt(v.taudif*p.deltat)*(2*DTE2[n] - v.heatflux_interior[n])
    
            v.heat_mixed[n] = v.heat_mixed[n-1] +v.heatflux_mixed[n] *(v.powtoheat*p.deltat)
    
            v.heat_interior[n] = v.heat_interior[n-1] + v.heatflux_interior[n] * (fso*v.powtoheat*p.deltat)
    
            v.temp[n] = flnd*v.temp_landair[n] + (1-flnd)*bsi*v.temp_sst[n]
        else
    
            # handle initial values
            DelQL = 0.0
            DelQ0 = 0.0
    
            QC1 = 0.0
            QC2 = 0.0
    
            DQ1 = 0.0
            DQ2 = 0.0
    
            DPAST2 = 0.0
            DPAST1 = 0.0
    
            v.temp_landair[TimestepIndex(1)] = 0.0
            v.temp_sst[TimestepIndex(1)] = 0.0
    
            v.heatflux_mixed[TimestepIndex(1)] = 0.0
            v.heatflux_interior[TimestepIndex(1)] = 0.0
    
            v.heat_mixed[TimestepIndex(1)] = 0.0
            v.heat_interior[TimestepIndex(1)] = 0.0
            v.temp[TimestepIndex(1)] = 0.0
        end
    end
    
end
