using MimiSNEASY
using Test

@testset "MimiSNEASY" begin

# TODO Enable these once they have been ported to Julia 1.0
# include("test_ccm.jl")
# include("test_doeclim.jl")
# include("test_sneasy.jl")

m = MimiSNEASY.get_model()
run(m)

end
