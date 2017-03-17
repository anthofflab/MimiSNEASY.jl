using DataFrames
using Gadfly
include("sneasy.jl")

# load data
df_emissions = readtable(joinpath(dirname(@__FILE__), "..", "data", "RCP85_EMISSIONS.csv"))
rename!(df_emissions, :YEARS, :year)
df_forcings = readtable(joinpath(dirname(@__FILE__), "..", "data", "forcing_rcp85.txt"), separator=' ')
df = join(df_emissions, df_forcings, on=:year)
df = DataFrame(year=df[:year], co2=df[:FossilCO2]+df[:OtherCO2], rf_aerosols=df[:aerosol_direct]+df[:aerosol_indirect], rf_other=df[:ghg_nonco2]+df[:solar]+df[:volcanic]+df[:other])

# Run SNEASY
init_fortran_sneasy()
MOC_strength, radiative_forc, ATM_CO2, atm_oc_flux, GL_surface_temp, GL_ocean_heat = run_fortran_sneasy(convert(Array, df[:co2]), convert(Array, df[:rf_other]), convert(Array, df[:rf_aerosols]), 2., 1.1, 1., 1.311, 0.502, 17.7, 0.047)
fin_fortran_sneasy()

# Add results to DataFrame
df[:MOC_strength] = MOC_strength
df[:radiative_forc] = radiative_forc
df[:ATM_CO2] = ATM_CO2
df[:atm_oc_flux] = atm_oc_flux
df[:GL_surface_temp] = GL_surface_temp
df[:GL_ocean_heat] = GL_ocean_heat

# Create empty directory for output
output_path = joinpath(dirname(@__FILE__), "..", "output", "demo_sneasy_jl")
if isdir(output_path)
	rm(output_path, recursive=true)
end
mkpath(output_path)

# Plot results
draw(PDF(joinpath(output_path, "MOC_strength.pdf"), 12cm, 12cm), plot(df, x="year", y="MOC_strength", Geom.line))
draw(PDF(joinpath(output_path, "radiative_forc.pdf"), 12cm, 12cm), plot(df, x="year", y="radiative_forc", Geom.line))
draw(PDF(joinpath(output_path, "ATM_CO2.pdf"), 12cm, 12cm), plot(df, x="year", y="ATM_CO2", Geom.line))
draw(PDF(joinpath(output_path, "atm_oc_flux.pdf"), 12cm, 12cm), plot(df, x="year", y="atm_oc_flux", Geom.line))
draw(PDF(joinpath(output_path, "GL_surface_temp.pdf"), 12cm, 12cm), plot(df, x="year", y="GL_surface_temp", Geom.line))
draw(PDF(joinpath(output_path, "GL_ocean_heat.pdf"), 12cm, 12cm), plot(df, x="year", y="GL_ocean_heat", Geom.line))
