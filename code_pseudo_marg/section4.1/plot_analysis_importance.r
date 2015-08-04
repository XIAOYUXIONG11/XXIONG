## Code to produce figure 2

pdf.options(width=10,height=10,pointsize=32)
ps.options(width=10,height=10,paper="special",horizontal=F, pointsize=32)

curvecolors = colors()[c(258, 26)]

## Choose the approximation and the number of data - RUN analysis_importance.r FIRST!!
APPROXIMATION = "EP"

for(N in c(50, 200))
{
## Curve #1
NPSEUDO = 1
filesave = paste("RES_PLOT_", "PROBIT", "_n_", N, "_method_", APPROXIMATION, "_nimp_", NPSEUDO, ".Rdata", sep=""); load(filesave)
c1.lq = exp.lq.mat.pseudo
c1.hq = exp.hq.mat.pseudo
c1.md = exp.median.mat.pseudo

## Curve #2
NPSEUDO = 64
filesave = paste("RES_PLOT_", "PROBIT", "_n_", N, "_method_", APPROXIMATION, "_nimp_", NPSEUDO, ".Rdata", sep=""); load(filesave)
c2.lq = exp.lq.mat.pseudo
c2.hq = exp.hq.mat.pseudo
c2.md = exp.median.mat.pseudo

## Generate the plot
fileplot = paste("PLOT_", "PROBIT", "_n_", N, "_method_", APPROXIMATION, ".eps", sep="")

if(APPROXIMATION == "LA") name.APPROXIMATION = "LA"
if(APPROXIMATION == "EP") name.APPROXIMATION = "EP"

postscript(fileplot)
par("mar"=c(2.9,2.9,1.1,0.4), "las"=0, "mgp"=c(1.8,0.6,0))
## Plot the median of the results for Nimp = 64
plot(exp(psi.tau.vect), c2.md, type="l", lwd=2, ylim = c(0,2), ylab="pseudo marginal", xlab="length-scale", col=1, main=paste("n = ", N, ", ", name.APPROXIMATION, sep=""))
## plot(exp(psi.tau.vect), c2.md, type="l", lwd=2, ylim = range(c(c1.lq, c2.lq, c1.hq, c2.hq)), ylab="pseudo marginal", xlab="length-scale", col=1, main=paste("n = ", N, ", ", name.APPROXIMATION, sep=""))

## Plot the prior
exp.psi.tau.vect = seq(0, max.tau.vect, by=0.01)
points(exp.psi.tau.vect, dgamma(exp.psi.tau.vect, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE), type="l", col=2, lwd=2) # "grey30")

## Curve 1
points(exp(psi.tau.vect), c1.lq, type="l", lwd=2, lty=5, col=curvecolors[1])
points(exp(psi.tau.vect), c1.hq, type="l", lwd=2, lty=5, col=curvecolors[1])

## Curve 2
points(exp(psi.tau.vect), c2.lq, type="l", lwd=2, lty=5, col=curvecolors[2])
points(exp(psi.tau.vect), c2.hq, type="l", lwd=2, lty=5, col=curvecolors[2])

## Legend
if(APPROXIMATION == "LA") {
  if(N == 50) legend(1.3, 1.8, legend=c("Nimp = 1", "Nimp = 64"), lwd=2, lty=5, col=c(curvecolors[1],curvecolors[2]))
  if(N == 200) legend(1.3, 1.8, legend=c("Nimp = 1", "Nimp = 64"), lwd=2, lty=5, col=c(curvecolors[1],curvecolors[2]))
}

if(APPROXIMATION == "EP") {
  if(N == 50) legend(1.3, 1.8, legend=c("Nimp = 1", "Nimp = 64"), lwd=2, lty=5, col=c(curvecolors[1],curvecolors[2]))
  if(N == 200) legend(1.3, 1.8, legend=c("Nimp = 1", "Nimp = 64"), lwd=2, lty=5, col=c(curvecolors[1],curvecolors[2]))
}

## True value of tau
abline(v=exp(TRUE.PSI.TAU), lty=2)
dev.off()
}
