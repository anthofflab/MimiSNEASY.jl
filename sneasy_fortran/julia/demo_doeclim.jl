using DataFrames
using Gadfly
include("sneasy.jl")

# load data
df = readtable(joinpath(dirname(@__FILE__), "..", "data", "forcing_rcp85.txt"), separator=' ')
df = DataFrame(year=df[:year], rf=df[:co2]+df[:aerosol_direct]+df[:aerosol_indirect]+df[:ghg_nonco2]+df[:solar]+df[:volcanic]+df[:other])

# Run doeclim
mod_temp, mod_heatflux_mixed, mod_heatflux_interior = run_fortran_doeclim(2., 1.1, convert(Array, df[:rf]))

# Add results to DataFrame
df[:temp] = mod_temp
df[:heatflux_mixed] = mod_heatflux_mixed
df[:heatflux_interior] = mod_heatflux_interior

# Create empty directory for output
output_path = joinpath(dirname(@__FILE__), "..", "output", "demo_doeclim_jl")
if isdir(output_path)
	rm(output_path, recursive=true)
end
mkpath(output_path)

# Plot results
draw(PDF(joinpath(output_path, "temp.pdf"), 12cm, 12cm), plot(df, x="year", y="temp", Geom.line))
draw(PDF(joinpath(output_path, "heatflux_mixed.pdf"), 12cm, 12cm), plot(df, x="year", y="heatflux_mixed", Geom.line))
draw(PDF(joinpath(output_path, "heatflux_interior.pdf"), 12cm, 12cm), plot(df, x="year", y="heatflux_interior", Geom.line))
