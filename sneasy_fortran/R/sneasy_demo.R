# demo of calling SNEASY Fortran model from R
# run R from command line (Unix command 'R'; may require 'module add R' on hammer)
# run this from R console:  source("sneasy_demo.R")

# R/Fortran interface is called in the sneasy() function (in file sneasy.R)

# before running, compile the Fortran source into a shared library
#  make sneasy.so


# load non-CO2 forcing time series.
forc_other = read.table("../sneasy/non_CO2_forcing.txt", col.names=c("year","Foth"))
F.oth = forc_other$Foth

# define and initialize DOECLIM model
source("sneasy.R")
is.loaded("../sneasy/sneasy")

# run SNEASY
setup.sneasy()

model.output = sneasy(S=2.0, kappa=1.1, alpha=0.6, Q10=1.311 , beta=0.502 , eta=17.722, hydsens=0.047, init.MOC = 22.6, init.CO2 = 285.2, endyear=mod.time[length(mod.time)])


# plot output
par(mfrow=c(3,2), mar=c(5,4,0.5,2))

plot(model.output$time, model.output$moc, col="blue", type="l", lwd=2, xlab="Year", ylab="MOC Strength [Sv]")

plot(model.output$time, model.output$forcing, col="red", type="l", lwd=2, xlab="Year", ylab="Radiative Forcing [W/m^2]")

plot(model.output$time, model.output$co2, col="green", type="l", lwd=2, xlab="Year", ylab="Atmospheric CO2 [ppm]")

plot(model.output$time, model.output$ocflux, col="red", type="l", lwd=2, xlab="Year", ylab="Atmosphere-Ocean flux [W/m^2]")

plot(model.output$time, model.output$temp, col="dark red", type="l", lwd=2, xlab="Year", ylab="Global Surface Temp. [K]")

plot(model.output$time, model.output$ocheat, col="magenta", type="l", lwd=2, xlab="Year", ylab="Global Ocean Heat [10^22 J]")


cleanup.sneasy()

