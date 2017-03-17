nsamp = 1000

sidx = sample(nrow(mcmc.out$batch), nsamp)
schain = mcmc.out$batch[sidx,]

pred.temp = array(dim=c(nsamp,length(mod.time)))
pred.ocheat = array(dim=c(nsamp,length(mod.time)))
pred.co2 = array(dim=c(nsamp,length(mod.time)))
pred.ocflux = array(dim=c(nsamp,length(mod.time)))
pred.moc = array(dim=c(nsamp,length(mod.time)))

for(i in 1:nsamp) {
	p = schain[i,]
	S = p[1]; kappa = p[2]; alpha = p[3]
	Q10 = p[4]; beta = p[5]; eta = p[6]
	hydsens = p[7]
	T0 = p[8]; H0 = p[9]; CO20 = p[10]; MOC0 = p[11]
	sigma.temp = p[12]; sigma.ocheat = p[13]; sigma.co2inst = p[14]; sigma.co2ice = p[15]
	rho.temp = p[16]; rho.ocheat = p[17]; rho.co2inst = p[18]
	
	model.output = sneasy(S, kappa, Q10, beta, eta, hydsens, endyear=2500)
	
	pred.temp[i,] = model.output$global.surf.temp + T0
	pred.ocheat[i,] = model.output$global.oc.heat + H0
	pred.co2[i,] = model.output$atmos.CO2 - model.output$atmos.CO2[1] + CO20
	pred.ocflux[i,] = model.output$atmos.oc.flux
	pred.moc[i,] = model.output$MOC.strength - model.output$MOC.strength[1] + MOC0	
}