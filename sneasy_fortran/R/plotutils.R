# plotutils.R
#
# Nathan M. Urban (nurban@psu.edu)
# Department of Geosciences, Penn State
#
# Miscellaneous plotting utilities.
#
# Plots Markov chains, marginal and pairwise joint distributions.
# The distribution plots can be histograms/scatterplots or
# smoothed 1D/2D kernel density estimates.  The pairwise
# plots also compute pairwise correlation coefficients.
#
# In the marginal distributions the red line is the mean and
# the blue dashed and dotted lines are the median and 2.5/97.5%
# quantiles (95% credible interval).  In the pair distributions
# the red line is a loess (locally weighted polynomial regression)
# smooth of the scatter (or contour) plot.


library(KernSmooth) # for kernel density estimation (smoothed histograms)

# open a new figure window for plotting
figure <- function(width=7, height=7)
{	
	if(! exists("nzchar"))
		nzchar = function(s){nchar(s)>0}
			
	defdev <- Sys.getenv("R_DEFAULT_DEVICE")
	if(!nzchar(defdev)) defdev <- "pdf"
	device <- if(interactive()) {
		intdev = Sys.getenv("R_INTERACTIVE_DEVICE")
		if(nzchar(intdev)) intdev
		else {
			if(.Platform$OS.type == "windows") "windows"
			else if (.Platform$GUI == "AQUA") "quartz"
			else if (Sys.getenv("DISPLAY") != "")
				switch(.Platform$GUI, "Tk" = "X11", "X11" = "X11", "GNOME" = "X11", defdev)
			else if (.Call("makeQuartzDefault")) "quartz"
			else defdev
		}
	} else defdev
    
	eval(parse(text=sprintf("%s(width=%d,height=%d)",device,width,height)))
}

# plot Markov chain
plot.chain <- function(chain, nthin=NA)
{
	if( !is.na(nthin) ) { # thin chain to have nthin elements
		chain <- chain[seq(1,dim(chain)[1],len=min(nthin,dim(chain)[1])),]
	}
	np <- dim(chain)[2]

	ncol <- ceiling(np/min(np,6)) # max 6 rows
	nrow <- ceiling(np/ncol)
	
	par(mfrow=c(nrow,ncol), ps=18, mar=c(5,5,4,2)+0.1, oma=c(1.5,0,0,0))
	#par(mfrow=c(nrow,ncol), ps=18, mar=c(5,5,4,2)+0.1, oma=c(0,0,1.2,0))

	for(i in 1:np) {
		par(mar=c(2,4.5,2,2))
		plot(chain[,i],type="l",ylab=colnames(chain)[i],frame=F)
	}
	
	#mtext("Markov chain", outer=TRUE, line=-1)	
}

plot.marginals <- function(chain, qlo=0.025, qhi=0.975, smoothing=1, truncate.kde=TRUE)
{
	N <- dim(chain)[1]; np <- dim(chain)[2]
	
	ncol <- floor(sqrt(np))
	nrow <- ceiling(np/ncol)
	par(mfrow=c(nrow,ncol))
	
	for(i in 1:np) {
		x <- chain[,i]
		est <- if(truncate.kde){bkde(x,bandwidth=smoothing*dpik(x),range.x=range(x))} else {bkde(x,bandwidth=smoothing*dpik(x))} # from KernSmooth
		#est <- if(truncate.kde, {density(x,adjust=smoothing,cut=0)} else {density(x,adjust=smoothing))}
		y <- est$y
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5*max(y)) )
		
		# plot pdf (smoothed kernel density estimate)
		plot.default(NA, xlim=range(x), ylim=c(0,max(y)), main=colnames(chain)[i], xlab=colnames(chain)[i], ylab="pdf", frame=F, cex.axis=1.2, cex.lab=1.2, cex.main=1.4)
		if(truncate.kde) {
			polygon(c(min(x),est$x,max(x)),c(0,y,0),col="gray")
		} else {
			polygon(est,col="gray")	
		}
		
		# plot mean and qlo/50/qhi quantiles
		lines(c(mean(x),mean(x)), c(0,1.05*max(y)), col="red")
		lines(c(median(x),median(x)), c(0,1.05*max(y)), col="blue", lty="dashed")
		lines(c(quantile(x,qlo),quantile(x,qlo)), c(0,1.05*max(y)), col="blue", lty="dotted")
		lines(c(quantile(x,qhi),quantile(x,qhi)), c(0,1.05*max(y)), col="blue", lty="dotted")	
	}
}

plot.pairs <- function(x, smooth=FALSE, smooth.hist=FALSE, smooth.scatter=FALSE, smoothing=1, truncate.kde=TRUE, qlo=0.025, qhi=0.975, bold=TRUE, ...)
{
	if(smooth) {
		smooth.hist = TRUE
		smooth.scatter = TRUE	
	}
	
	histpanel <- function(x,smooth,...) {panel.hist(x,smooth.hist,smoothing,truncate.kde,qlo,qhi,...)}
	scatterpanel <- function(x,y,smooth,...) {panel.scatter(x,y,smooth.scatter,smoothing,truncate.kde,...)}
	
	pairs(x, lower.panel=scatterpanel, diag.panel=histpanel, upper.panel=panel.cor, font.labels=ifelse(bold,2,1))
}

panel.hist <- function(x, smooth, smoothing=1, truncate.kde=TRUE, qlo=0.025, qhi=0.975, ...)
{
	if(smooth) { # kernel density estimate
		est <- ifelse(truncate.kde, bkde(x,bandwidth=smoothing*dpik(x),range.x=range(x)), bkde(x,bandwidth=smoothing*dpik(x))) # from KernSmooth
		#est <- ifelse(truncate.kde, density(x,adjust=smoothing,cut=0), density(x,adjust=smoothing))
		y <- est$y
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5*max(y)) )
		if(truncate.kde) {
			polygon(c(min(x),est$x,max(x)),c(0,y,0),col="gray")
		} else {
			polygon(est,col="gray")	
		}
	} else { # histogram
		h <- hist(x, plot = FALSE)
		breaks <- h$breaks; nB <- length(breaks)
    	y <- h$density
    	usr <- par("usr"); on.exit(par(usr))
    	par(usr = c(usr[1:2], 0, 1.5*max(y)) )
    	rect(breaks[-nB], 0, breaks[-1], y, border=NA, col="gray", ...)
    	n <- nB-1
    	outline.x0 <- rep(NA, 2*n-1)
    	outline.y0 <- rep(NA, 2*n-1)
    	outline.x1 <- rep(NA, 2*n-1)
    	outline.y1 <- rep(NA, 2*n-1)
    	for(i in 1:n) {
    		outline.x0[2*i-1] <- breaks[i]
    		outline.y0[2*i-1] <- y[i]
    		outline.x1[2*i-1] <- breaks[i+1]
    		outline.y1[2*i-1] <- y[i]
    		if(i < n) {
    			outline.x0[2*i] <- breaks[i+1]
    			outline.y0[2*i] <- y[i]
    			outline.x1[2*i] <- breaks[i+1]
    			outline.y1[2*i] <- y[i+1]
    		}
    	}
    	segments(outline.x0,outline.y0,outline.x1,outline.y1)
	}
	
	# plot mean and qlo/50/qhi quantiles
	lines(c(mean(x),mean(x)), c(0,1.05*max(y)), col="red")
	lines(c(median(x),median(x)), c(0,1.05*max(y)), col="blue", lty="dashed")
	lines(c(quantile(x,qlo),quantile(x,qlo)), c(0,1.05*max(y)), col="blue", lty="dotted")
	lines(c(quantile(x,qhi),quantile(x,qhi)), c(0,1.05*max(y)), col="blue", lty="dotted")
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * r, col="blue")
}

panel.scatter <- function (x, y, smooth, smoothing=1, col = par("col"), bg = NA, pch = 20, cex = 0.5, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
	if(smooth) { # 2D kernel density estimate
		maxpts <- min(50000,length(x))
		idx <- seq(1,maxpts,length=maxpts); x <- x[idx]; y <- y[idx]
		z <- cbind(x,y)
		est <- bkde2D(z, bandwidth=smoothing*c(dpik(x),dpik(y))) # from KernSmooth
		l <- contourLines(est$x1, est$x2, est$fhat, levels=seq(0.05,0.95,0.15)*max(est$fhat))
		for(i in 1:length(l)) {
			lines(l[[i]]$x,l[[i]]$y)
		}
		ok <- is.finite(x) & is.finite(y)
    	if (any(ok)) 
       	lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
	} else { # scatterplot
		# subsample points to prevent scatterplot from being too crowded
		nsamp <- 1000
		maxpts <- min(nsamp, length(x))
		idx <- sample(1:max(length(x),maxpts), maxpts); xs <- x[idx]; ys <- y[idx]
		points(xs, ys, pch = pch, col = "black", bg = bg, cex = cex)
		# ... but compute smoothed curve using full data set
    	ok <- is.finite(xs) & is.finite(ys)
    	if (any(ok)) 
        	lines(stats::lowess(xs[ok], ys[ok], f = span, iter = iter), 
            	col = col.smooth, ...)	
	}
}