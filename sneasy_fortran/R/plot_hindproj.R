par(mfrow=c(3,2), mar=c(3.5,3.5,1.5,0), mgp=c(2.3,0.9,0), cex.lab=1.3, cex.axis=1.2, cex.main=1.3)

plot(NA, xlim=range(mod.time), ylim=c(-3,12), xlab="Year", ylab="Forcing", main="Radiative forcing", frame=FALSE)
polygon(c(mod.time,rev(mod.time)), c(q.rf.lo,rev(q.rf.hi)), col="gray", border=NA)
lines(mod.time, mean.rf, lwd=2, col="purple")

plot(NA, xlim=range(mod.time), ylim=c(280,2000), xlab="Year", ylab=expression(CO[2]), main=expression(bold(paste("Atmospheric ",CO[2]," concentration",sep=""))), frame=FALSE)
polygon(c(mod.time,rev(mod.time)), c(q.co2.lo,rev(q.co2.hi)), col="gray", border=NA)
lines(mod.time, mean.co2, lwd=2, col="purple")
points(obs.co2ice.time, obs.co2ice, pch=19, cex=0.5, col="red")
points(obs.co2inst.time, obs.co2inst, pch=19, cex=0.5, col="red")

plot(NA, xlim=range(mod.time), ylim=c(-0.5,10), xlab="Year", ylab="Temperature", main="Surface temperature anomaly", frame=FALSE)
polygon(c(mod.time,rev(mod.time)), c(q.temp.lo,rev(q.temp.hi)), col="gray", border=NA)
lines(mod.time, mean.temp, lwd=2, col="purple")
points(obs.temp.time, obs.temp, pch=19, cex=0.5, col="red")

plot(NA, xlim=range(mod.time), ylim=c(-100,2500), xlab="Year", ylab="Ocean heat", main="Ocean heat anomaly", frame=FALSE)
polygon(c(mod.time,rev(mod.time)), c(q.ocheat.lo,rev(q.ocheat.hi)), col="gray", border=NA)
lines(mod.time, mean.ocheat, lwd=2, col="purple")
points(obs.ocheat.time, obs.ocheat, pch=19, cex=0.5, col="red")

plot(NA, xlim=range(mod.time), ylim=c(-10,5), xlab="Year", ylab="Ocean C flux", main="Atmosphere-ocean carbon flux", frame=FALSE)
polygon(c(mod.time,rev(mod.time)), c(q.ocflux.lo,rev(q.ocflux.hi)), col="gray", border=NA)
lines(mod.time, mean.ocflux, lwd=2, col="purple")
points(obs.ocflux.time, obs.ocflux, pch=19, cex=0.5, col="red")

plot(NA, xlim=range(mod.time), ylim=c(0,25), xlab="Year", ylab="MOC", main="MOC strength", frame=FALSE)
polygon(c(mod.time,rev(mod.time)), c(q.moc.lo,rev(q.moc.hi)), col="gray", border=NA)
lines(mod.time, mean.moc, lwd=2, col="purple")
points(obs.moc.time, obs.moc, pch=19, cex=0.5, col="red")