# Updated to call 'fBasics' instead of 'fUtilities'
# Last modified Fri 18 May 2012 @ 10:52:29 EDT by Rob Nicholas.


# hindcasts and projections
source("loadobs.R")
source("sneasy.R")

# simulate stationary AR(1) process
ar1.sim = function(N,rho1,sigma) {
	x = rep(NA,N)
	x[1] = sigma/sqrt(1-rho1^2)
	for(i in 2:N)
		x[i] = rho1*x[i-1] + rnorm(1,sd=sigma)
	return(x)
}

# calculate total forcing for a given CO2 concentration and aerosol forcing scale factor
total.forcing = function(co2, alpha) 3.7*log(co2/co2[1])/log(2) + nonco2.forcing(forcing, alpha)

# thinned (short) chain
nsamp = 10000
sidx = seq(1,nrow(chain),len=nsamp) # thin chain further to 10000 samples
schain = chain[sidx,]

N = length(mod.time)
rf = array(dim=c(nrow(schain),N))
pred.temp.mod = array(dim=c(nrow(schain),N))
pred.ocheat.mod = array(dim=c(nrow(schain),N))
pred.co2.mod = array(dim=c(nrow(schain),N))
pred.ocflux.mod = array(dim=c(nrow(schain),N))
pred.moc.mod = array(dim=c(nrow(schain),N))
pred.temp.noise = array(dim=c(nrow(schain),N))
pred.ocheat.noise = array(dim=c(nrow(schain),N))
pred.co2.noise = array(dim=c(nrow(schain),N))
pred.ocflux.noise = array(dim=c(nrow(schain),N))
pred.moc.noise = array(dim=c(nrow(schain),N))

setup.sneasy()

mod.best = sneasy(p0[1],p0[2],p0[3],p0[4],p0[5],p0[6],p0[7],p0[10],p0[11])
temp.best = mod.best$temp + p0[8]
ocheat.best = mod.best$ocheat + p0[9]
co2.best = mod.best$co2
ocflux.best = mod.best$ocflux
moc.best = mod.best$moc

# generate hindcasts 
for(i in 1:nrow(schain)) {
	mod.out = sneasy(schain[i,1],schain[i,2],schain[i,3],schain[i,4],schain[i,5],schain[i,6],schain[i,7],schain[i,10],schain[i,11])

	pred.temp.mod[i,] = mod.out$temp + schain[i,8]
	pred.ocheat.mod[i,] = mod.out$ocheat + schain[i,9]
	pred.co2.mod[i,] = mod.out$co2
	pred.ocflux.mod[i,] = mod.out$ocflux
	pred.moc.mod[i,] = mod.out$moc
	
	rf[i,] = total.forcing(pred.co2.mod[i,], schain[i,3])

	pred.temp.noise[i,] = pred.temp.mod[i,] + ar1.sim(N,schain[i,16],schain[i,12])
	pred.ocheat.noise[i,] = pred.ocheat.mod[i,] + ar1.sim(N,schain[i,17],schain[i,13])
	pred.co2.noise[i,] = pred.co2.mod[i,] + ar1.sim(N,schain[i,18],schain[i,14])
	pred.ocflux.noise[i,] = pred.ocflux.mod[i,] + rnorm(N,0,0.4*sqrt(10)) # convert decadal to annual error
	pred.moc.noise[i,] = pred.moc.mod[i,] + rnorm(N,0,2.0) # +/- 2 Sv
}

cleanup.sneasy()

# mean, median, ci hindcasts
# The 'fUtilities' library has been withdrawn from CRAN; use 'fBasics' instead
library(fBasics)

mean.rf = colMeans(rf)
sd.rf = colSds(rf)
skew.rf = colSkewness(rf)
kurt.rf = colKurtosis(rf)
q.rf.lo = colQuantiles(rf, 0.025)
q.rf.med = colQuantiles(rf, 0.5)
q.rf.hi = colQuantiles(rf, 0.975)

mean.temp = colMeans(pred.temp.noise)
sd.temp = colSds(pred.temp.noise)
skew.temp = colSkewness(pred.temp.noise)
kurt.temp = colKurtosis(pred.temp.noise)
q.temp.lo = colQuantiles(pred.temp.noise, 0.025)
q.temp.med = colQuantiles(pred.temp.noise, 0.5)
q.temp.hi = colQuantiles(pred.temp.noise, 0.975)

mean.ocheat = colMeans(pred.ocheat.noise)
sd.ocheat = colSds(pred.ocheat.noise)
skew.ocheat = colSkewness(pred.ocheat.noise)
kurt.ocheat = colKurtosis(pred.ocheat.noise)
q.ocheat.lo = colQuantiles(pred.ocheat.noise, 0.025)
q.ocheat.med = colQuantiles(pred.ocheat.noise, 0.5)
q.ocheat.hi = colQuantiles(pred.ocheat.noise, 0.975)

mean.co2 = colMeans(pred.co2.noise)
sd.co2 = colSds(pred.co2.noise)
skew.co2 = colSkewness(pred.co2.noise)
kurt.co2 = colKurtosis(pred.co2.noise)
q.co2.lo = colQuantiles(pred.co2.noise, 0.025)
q.co2.med = colQuantiles(pred.co2.noise, 0.5)
q.co2.hi = colQuantiles(pred.co2.noise, 0.975)

mean.ocflux = colMeans(pred.ocflux.noise)
sd.ocflux = colSds(pred.ocflux.noise)
skew.ocflux = colSkewness(pred.ocflux.noise)
kurt.ocflux = colKurtosis(pred.ocflux.noise)
q.ocflux.lo = colQuantiles(pred.ocflux.noise, 0.025)
q.ocflux.med = colQuantiles(pred.ocflux.noise, 0.5)
q.ocflux.hi = colQuantiles(pred.ocflux.noise, 0.975)

mean.moc = colMeans(pred.moc.noise)
sd.moc = colSds(pred.moc.noise)
skew.moc = colSkewness(pred.moc.noise)
kurt.moc = colKurtosis(pred.moc.noise)
q.moc.lo = colQuantiles(pred.moc.noise, 0.025)
q.moc.med = colQuantiles(pred.moc.noise, 0.5)
q.moc.hi = colQuantiles(pred.moc.noise, 0.975)

if(save.output) {
	fstr = "_hindproj"
	save(nsamp, sidx, schain, rf, pred.temp.mod, pred.ocheat.mod, pred.co2.mod, pred.ocflux.mod, pred.moc.mod, pred.temp.noise, pred.ocheat.noise, pred.co2.noise, pred.ocflux.noise, pred.moc.noise, file=paste("../output/pred",fstr,"_o",atype,"_y",endyear,".Rdata",sep=""))
	save(nsamp, sidx, schain, mod.time, mean.rf, sd.rf, skew.rf, kurt.rf, q.rf.lo, q.rf.med, q.rf.hi, mean.temp, sd.temp, skew.temp, kurt.temp, q.temp.lo, q.temp.med, q.temp.hi, mean.ocheat, sd.ocheat, skew.ocheat, kurt.ocheat, q.ocheat.lo, q.ocheat.med, q.ocheat.hi, mean.co2, sd.co2, skew.co2, kurt.co2, q.co2.lo, q.co2.med, q.co2.hi, mean.ocflux, sd.ocflux, skew.ocflux, kurt.ocflux, q.ocflux.lo, q.ocflux.med, q.ocflux.hi, mean.moc, sd.moc, skew.moc, kurt.moc, q.moc.lo, q.moc.med, q.moc.hi, file=paste("../output/statpred",fstr,"_o",atype,"_y",endyear,".Rdata",sep=""))
}