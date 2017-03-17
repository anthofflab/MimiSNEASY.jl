using DataFrames
using Gadfly
include("sneasy.jl")

# load data
df = readtable(joinpath(dirname(@__FILE__), "..", "data", "RCP85_EMISSIONS.csv"))
rename!(df, :YEARS, :year)
df = DataFrame(year=df[:year], co2=df[:FossilCO2]+df[:OtherCO2])

# Run ccm
init_fortran_ccm()
atmco2_out = run_fortran_ccm(2., 1., 1.311, 0.502, 280., ones(736).*2.1, convert(Array, df[:co2]))
fin_fortran_ccm()

# Add results to DataFrame
df[:atmco2_out] = atmco2_out

# Create empty directory for output
output_path = joinpath(dirname(@__FILE__), "..", "output", "demo_ccm_jl")
if isdir(output_path)
	rm(output_path, recursive=true)
end
mkpath(output_path)

# Plot results
draw(PDF(joinpath(output_path, "atmco2_out.pdf"), 12cm, 12cm), plot(df, x="year", y="atmco2_out", Geom.line))
