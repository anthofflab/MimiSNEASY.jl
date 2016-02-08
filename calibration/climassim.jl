using Lora
using Gadfly
using NLopt

# Should the mcmc algorithm start with the same start values as the R code?
strict_R_compat = true
# Multiple chains? This was in the R version but doesn't currently work here
multiple_chains = false
# Possible falues are :fortran and :julia
sneasy_version = :julia

include("assim.jl")

run_mimi_sneasy! = construct_run_mimi_sneasy()

# Get log-likelihood function
if sneasy_version==:fortran
    log_post = construct_log_post(run_fortran_sneasy!)
elseif sneasy_version==:julia
    log_post = construct_log_post(run_mimi_sneasy!)
else
    error("Unknown sneasy version requested")
end

parnames = ["S","kappa","alpha","Q10","beta","eta","h","T0","H0","CO20","MOC0","sigma.T","sigma.H","sigma.CO2.inst","sigma.CO2.ice","rho.T","rho.H","rho.CO2.inst"]

initial_guess = [2.7,2.9,1.0,4.2,0.9,23,0.03,-0.06,-33,286,19.5,0.1,2.0,0.45,2.25,0.55,0.9,0.95]

# This is for Frank to get working first
# log_post(initial_guess)

p0 = zeros(18)
if strict_R_compat
	# These are the parameter estimates that R gets from maximising the posterior
	p0 = [2.67391818, 3.13600522,   0.98236787,   4.17150643,   0.88898206,  23.01048615,   0.04972258,  -0.03610151, -25.70232473, 285.91278519,  19.71299613,   0.09730246, 2.04538221,   0.44754775,   2.73731543,   0.49014552,   0.81813141,   0.97858006]
else
	opt = Opt(:LN_NELDERMEAD, 18)
	lower_bounds!(opt, Float64[0,0,0,0,0,0,0,-Inf,-100,280,10,0,0,0,0,0,0,0])
	upper_bounds!(opt, Float64[Inf,Inf,3,5,1,200,0.06,Inf,0,295,30,0.2,4,1,10,0.99,0.99,0.99])
	max_objective!(opt, (x, grad)->log_post(x))
	maxtime!(opt, 180)

	(minf,minx,ret) = optimize(opt, initial_guess)
	p0 = minx
end

model = likelihood_model(BasicContMuvParameter(:p, logtarget=log_post), false)

step = [1.6,1.7,0.25,0.75,0.15,40,0.015,0.03,9,0.7,1.3,0.005,0.25,0.045,0.57,0.07,0.06,0.11]./10

if multiple_chains
    job1 = BasicMCJob(model, MH(step), BasicMCRange(nsteps=10000, burnin=0), Dict(:p=>p0))
    @time mcchain1 = run(job1)

    job2 = BasicMCJob(model, MH(proposal_matrix(output(job1),mult=0.5)), BasicMCRange(nsteps=100_000, burnin=0), Dict(:p=>p0))
    @time mcchain2 = run(job2)

    job3 = BasicMCJob(model, MH(proposal_matrix(output(job2),mult=0.5)), BasicMCRange(nsteps=1_000_000, burnin=0), Dict(:p=>p0))
    @time mcchain3 = run(job3)

    job = job3
else
    job = BasicMCJob(model, MH(step), BasicMCRange(nsteps=100_000, burnin=1_000), Dict(:p=>p0))
    @time mcchain1 = run(job)
end

for (i,name) in enumerate(parnames)
    draw(PDF("plot-$name.pdf", 12cm,12cm), plot(x=vec(output(job).value[i,:]), Geom.density, Guide.title(name)))
end

#cleanup_sneasy() # deallocates memory after SNEASY is done

#println("Acceptance: $(acceptance(mcchain3))") # acceptance rates around 0.25 (+/- 0.05) are good
