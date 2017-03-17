# Updated to call 'fBasics' instead of 'fUtilities'
# Last modified Fri 18 May 2012 @ 10:50:27 EDT by Rob Nicholas.

save.output = TRUE

assim.temp = TRUE
assim.ocheat = TRUE
assim.co2inst = TRUE
assim.co2ice = TRUE
assim.ocflux = FALSE
assim.moc = FALSE

endyear = 2010

source("sneasy.R") # SNEASY Earth system model
source("loadobs.R") # load observations
source("assim.R") # Bayesian parameter estimation

setup.sneasy() # call this to initialize SNEASY

# parameters: climate sensitivity, ocean vertical diffusivity, aerosol forcing scale factor, initial surface temperature, initial ocean heat, surface temperature residual error, ocean heat residual error, surface temperature residual autocorrelation, ocean heat residual autocorrelation
parnames = c("S","kappa","alpha","Q10","beta","eta","h","T0","H0","CO20","MOC0","sigma.T","sigma.H","sigma.CO2.inst","sigma.CO2.ice","rho.T","rho.H","rho.CO2.inst")

# "best guess" parameter values to initialize Markov chain
# start with hard-coded guess, then maximize posterior with optim()
# (can use DEoptim() instead of optim() for global optimization)
#p0 = optim(c(3.7,3.6,1.1,4.2,0.9,23,0.03,-0.06,-33,286,19.5,0.1,2.0,0.45,2.25,0.55,0.9,0.95), function(p) -log.post(p))$par
p0 = c(2.7,2.9,1.0,4.2,0.9,23,0.03,-0.06,-33,286,19.5,0.1,2.0,0.45,2.25,0.55,0.9,0.95)
p0 = optim(p0, function(p) -log.post(p))$par
print(p0)
names(p0) = parnames


	##### assimilate observations
# (two shorter assimilations to adaptively estimate posterior
# covariance, then a longer final assimilation)

library(mcmc)

	t0 = Sys.time()

step = c(1.6,1.7,0.25,0.75,0.15,40,0.015,0.03,9,0.7,1.3,0.005,0.25,0.045,0.57,0.07,0.06,0.11)/10

mcmc.out1 = metrop(log.post, p0, 10000, scale=step)
prechain1 = mcmc.out1$batch

mcmc.out2 = metrop(log.post, p0, 100000, scale=proposal.matrix(prechain1,mult=0.5))
prechain2 = mcmc.out2$batch

mcmc.out = metrop(log.post, p0, 1000000, scale=proposal.matrix(prechain2,mult=0.5))

cleanup.sneasy() # deallocates memory after SNEASY is done

	print("MCMC run time:"); print(Sys.time()-t0)

chain = mcmc.out$batch
colnames(chain) = parnames

sprintf("MCMC acceptance rate = %.3f", mcmc.out$accept) # acceptance rates around 0.25 (+/- 0.05) are good
print("posterior summary:")
print(summary(chain))

p.best = p0 # if the MCMC algorithm kept track of the maximum posterior point in the chain, should use that instead of the start point
print("best fit parameters:")
print(p.best)

# compute statistics
# The 'fUtilities' library has been withdrawn from CRAN; use 'fBasics' instead
library(fBasics)

chain.range = apply(chain, 2, range)
chain.density = apply(chain, 2, function(x) density(x, adj=2))
chain.cdf = apply(chain, 2, ecdf)
chain.q = apply(chain, 2, function(x) quantile(x, c(0.001,seq(0.01,0.99,by=0.01),0.999)))

chain.mean = colMeans(chain)
chain.sd = colSds(chain)
chain.skew = colSkewness(chain)
chain.kurt = colKurtosis(chain)
chain.lo = colQuantiles(chain, 0.025)
chain.med = colQuantiles(chain, 0.5)
chain.hi = colQuantiles(chain, 0.975)

chain.SIPCC = sum(chain[,1]>2 & chain[,1]<4.5)/nrow(chain)
chain.S15 = sum(chain[,1]>1.5)/nrow(chain)
chain.S3 = sum(chain[,1]>3)/nrow(chain)
chain.S45 = sum(chain[,1]>4.5)/nrow(chain)
chain.S6 = sum(chain[,1]>6)/nrow(chain)

if(save.output) {
	atype = ""
	if(assim.temp) atype = paste(atype,"T",sep="")
	if(assim.ocheat) atype = paste(atype,"H",sep="")
	
	# thin chain to 100000 samples before saving
	chain = mcmc.out$batch[seq(1,nrow(mcmc.out$batch),len=100000),]
	save(chain, file=paste("../output/chain_hist_o",atype,"_y",endyear,".Rdata",sep=""))
	
	save(chain.range, chain.density, chain.cdf, chain.q, chain.mean, chain.sd, chain.skew, chain.kurt, chain.lo, chain.med, chain.hi, chain.S15, chain.S3, chain.S45, chain.S6, file=paste("../output/statchain_hist_o",atype,"_y",endyear,".Rdata",sep=""))
}

source("hindproj.R")

dev.new()
source("plot_hindproj.R")

dev.new()
source("plot_pdfs.R")