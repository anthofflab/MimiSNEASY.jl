# plot_pdfsMAP.R    Make a plot of pdfs as is done in climassim.R with
# MAP estimate for each indicated by a vertical line.

assim.temp = TRUE
assim.ocheat = TRUE
assim.co2inst = TRUE
assim.co2ice = TRUE
assim.ocflux = TRUE
assim.moc = TRUE

# Load workspace or just variables from files.
load("../output/statchain_hist_oTH_y2010.Rdata")
load("../output/chain_hist_oTH_y2010.Rdata")

endyear = 2010

parnames = c("S","kappa","alpha","Q10","beta","eta","h","T0","H0","CO20","MOC0","sigma.T","sigma.H","sigma.CO2.inst","sigma.CO2.ice","rho.T","rho.H","rho.CO2.inst")

source("sneasy.R")
source("assim.R")
source("loadobs.R")

setup.sneasy()

# Do the optimization
p.best = optim(chain.mean, function(p) -log.post(p))$par

cleanup.sneasy()

dev.new()
par(mfrow=c(5,4), mar=c(3,2.5,1.5,0.5), cex.lab=1.3, cex.axis=1.2, cex.main=1.4)

for(i in 1:ncol(chain)) {
	plot(density(chain[,i]), type="h", col="gray", lwd=2, xlab="", ylab="", main=parnames[i], frame=FALSE)
	lines(density(chain[,i]))
    abline(v=p.best[[i]])
}

pdf(file="../figures/marginals_MAP.pdf")
par(mfrow=c(5,4), mar=c(3,2.5,1.5,0.5), cex.lab=1.3, cex.axis=1.2, cex.main=1.4)

for(i in 1:ncol(chain)) {
	plot(density(chain[,i]), type="h", col="gray", lwd=2, xlab="", ylab="", main=parnames[i], frame=FALSE)
	lines(density(chain[,i]))
    abline(v=p.best[[i]])
}
dev.off()
