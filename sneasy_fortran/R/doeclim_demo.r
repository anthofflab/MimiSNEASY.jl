# demo of calling DOECLIM Fortran model from R
# run R from command line (Unix command 'R'; may require 'module add R' on hammer)
# run this from R console:  source("doeclim_demo.R")

# R/Fortran interface is called in the doeclim() function (in file doeclim.R)

# before running, compile the Fortran source into a shared library
#  option 1:  R CMD SHLIB doeclim.f90
#  option 2:  gfortran -dynamiclib doeclim.f90 -o doeclim.so


# define and initialize DOECLIM model
source("doeclim.R")


# run DOECLIM
model.output = doeclim(S=2, kappa=1.1, alpha=0.6)


# plot output
par(mfrow=c(2,1), mar=c(5,4,0.5,2))

plot(model.output$time, model.output$temp, col="red", type="l", lwd=2, xlab="Year", ylab="Temperature (K)")

plot(model.output$time, model.output$ocheat, col="blue", type="l", lwd=2, xlab="Year", ylab="Ocean heat (10^22 J)")