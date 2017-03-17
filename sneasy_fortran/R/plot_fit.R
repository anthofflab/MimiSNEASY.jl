p = optim(c(3.6,3.7,1.1,4.2,0.9,22,0.02,-0.05,-33,286,19,0.1,2.1,0.4,2.3,0.6,0.9,0.95), function(p) -log.post(p))$par

	S = p[1]; kappa = p[2]; alpha = p[3]
	Q10 = p[4]; beta = p[5]; eta = p[6]
	hydsens = p[7]
	T0 = p[8]; H0 = p[9]; CO20 = p[10]; MOC0 = p[11]
	sigma.temp = p[12]; sigma.ocheat = p[13]; sigma.co2inst = p[14]; sigma.co2ice = p[15]
	rho.temp = p[16]; rho.ocheat = p[17]; rho.co2inst = p[18]
	
model.out = sneasy(S, kappa, Q10, beta, eta, hydsens)

par(mfrow=c(3,2), mar=c(4,4,1,1))

plot(obs.temp.time, obs.temp, col="red")
lines(mod.time, model.out$global.surf.temp+T0, lwd=2)

plot(obs.ocheat.time, obs.ocheat, col="red", ylim=c(-15,25))
lines(mod.time, model.out$global.oc.heat+H0, lwd=2)

plot(obs.co2inst.time, obs.co2inst, col="red", xlim=c(1850,2010), ylim=c(280,390))
points(obs.co2ice.time, obs.co2ice, col="red")
lines(mod.time, model.out$atmos.CO2, lwd=2)

plot(obs.ocflux.time, obs.ocflux, col="red", xlim=c(1850,2010), ylim=c(-2.5,0))
lines(mod.time, model.out$atmos.oc.flux, lwd=2)

plot(obs.moc.time, obs.moc, col="red")
lines(mod.time, model.out$MOC.strength-model.out$MOC.strength[1]+MOC0, lwd=2)