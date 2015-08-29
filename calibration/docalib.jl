using Mimi
using Lora

include("assim.jl")

logli = construct_loglikelihood()

initial_guess = [2.7,2.9,1.0,4.2,0.9,23,0.03,-0.06,-33,286,19.5,0.1,2.0,0.45,2.25,0.55,0.9,0.95]

mcmodel = model(logli, init=initial_guess)

step = [1.6,1.7,0.25,0.75,0.15,40,0.015,0.03,9,0.7,1.3,0.005,0.25,0.045,0.57,0.07,0.06,0.11]./10

mcchain = Lora.run(mcmodel, MH(x::Vector{Float64} -> rand(MvNormal(x, step))), SerialMC(nsteps=100_000, burnin=5000))

