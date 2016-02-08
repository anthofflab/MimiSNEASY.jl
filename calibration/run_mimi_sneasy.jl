include("../src/sneasy.jl")

function construct_run_mimi_sneasy()
    m = getsneasy(nsteps=246)

    function run_mimi_sneasy!(
	       MOC_strength::Vector{Float64},
           radiative_forc::Vector{Float64},
           ATM_CO2::Vector{Float64},
           atm_oc_flux::Vector{Float64},
           GL_surface_temp::Vector{Float64},
           GL_ocean_heat::Vector{Float64},
           co2_emissions::Vector{Float64},
           rf_aerosol::Vector{Float64},
           rf_nonco2::Vector{Float64},
           S::Float64,
           kappa::Float64,
           alpha::Float64,
           Q10::Float64,
           beta::Float64,
           eta::Float64,
           hydsens::Float64,
           init_CO2::Float64=280.,
           init_MOC::Float64=20.)

        n = length(co2_emissions)

        if length(rf_aerosol)!=n || length(rf_nonco2)!=n
            error("All input vectors must have the same length")
        end

        setparameter(m, :doeclim, :t2co, S)
        setparameter(m, :doeclim, :kappa, kappa)
        setparameter(m, :ccm, :Q10, Q10)
        setparameter(m, :ccm, :Beta, beta)
        setparameter(m, :ccm, :Eta, eta)
        setparameter(m, :ccm, :atmco20, init_CO2)
        setparameter(m, :radforc, :alpha, alpha)

        setparameter(m, :ccm, :CO2_emissions, co2_emissions)
        setparameter(m, :radforc, :forcing_other, zeros(n))
        setparameter(m, :radforc, :forcing_volcanic, zeros(n))
        setparameter(m, :radforc, :forcing_solar, zeros(n))
        setparameter(m, :radforc, :forcing_ghg_nonco2, rf_nonco2)
        setparameter(m, :radforc, :forcing_aerosol_direct, rf_aerosol)
        setparameter(m, :radforc, :forcing_aerosol_indirect, zeros(n))

        run(m)

        MOC_strength[:] = NaN
        radiative_forc[:] = m[:radforc, :rf]
        ATM_CO2[:] = m[:ccm, :atmco2]
        atm_oc_flux[:] = m[:ccm, :atm_oc_flux]
        GL_surface_temp[:] = m[:doeclim, :temp]
        GL_ocean_heat[:] = m[:doeclim, :heat_interior]

        return
    end

    return run_mimi_sneasy!
end
