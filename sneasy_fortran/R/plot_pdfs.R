par(mfrow=c(5,4), mar=c(3,2.5,1.5,0.5), cex.lab=1.3, cex.axis=1.2, cex.main=1.4)

for(i in 1:ncol(chain)) {
	plot(density(chain[,i]), type="h", col="gray", lwd=2, xlab="", ylab="", main=parnames[i], frame=FALSE)
	lines(density(chain[,i]))
}