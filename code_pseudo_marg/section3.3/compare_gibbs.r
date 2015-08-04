## Code to generate the plot in figure 1
## Compare different parametrizations - Ancillary Augmentation (AA), Sufficient Augmentation (SA), and Surrogate method (SURR)

## One hyperparameter only

## **************************************************
## ******************** Initial parameters
## **************************************************

N = 200
NPSEUDO = 100

## APPROXIMATION = "LA"
APPROXIMATION = "LA"

COVARIANCE.MODE = "H"
NREP = 10
D = 2

## **************************************************
## ******************** Load useful functions
## **************************************************

## ps.options(width=18,height=18,paper="special",horizontal=F, pointsize=32)

pdf.options(width=10,height=10,pointsize=32)
ps.options(width=20,height=10,paper="special",horizontal=F, pointsize=32)

APPROACH = "PM"
source("../functions/PRIORS.r")
source("../functions/GP.FUNCTIONS.r")
source("../functions/EP.r")
source("../functions/LA.r")

reset.number.calls()

source("../functions/GENERATE.DATA.r")

MODEL = new.env()
MODEL$NLS = NLS
MODEL$SIZE.GROUPS = c(1, MODEL$NLS)
names(MODEL$SIZE.GROUPS) = MODEL$NAMES.GROUPS = c("psi.sigma", "psi.tau")

if(APPROXIMATION == "LA") compute.approximation = LA
if(APPROXIMATION == "EP") compute.approximation = EP.NOGRAD


## ***************************************************
## ******************** Loop through the hyperparameter and compute the pseudo marginal
## ***************************************************


set.seed(1)

CURRENT.STATE <<- new.env()

min.tau.vect = 0.01
if(N == 50) max.tau.vect = 2.8
if(N == 200) max.tau.vect = 2.8

NPSITAU = 80
psi.tau.vect = log(seq(min.tau.vect, max.tau.vect, length.out=NPSITAU))
step.psi.tau.vect = diff(seq(min.tau.vect, max.tau.vect, length.out=NPSITAU)[c(1,2)])
mat.pseudo = matrix(0, NREP, length(psi.tau.vect))


for(iii in 1:length(psi.tau.vect))
  {
    cat(formatC(iii/length(psi.tau.vect)*100, format="f", digits=2), "%    \r")
    CURRENT.STATE$psi.tau = psi.tau.vect[iii]
    CURRENT.STATE$psi.sigma = log(TRUE.SIGMA)

    compute.log.p.theta(CURRENT.STATE)

    compute.L.chol.Q.mat(CURRENT.STATE)

    for(jjj in 1:NREP)
      {
        compute.log.p.y.giv.theta.pseudo.marginal(CURRENT.STATE)
        
        
        mat.pseudo[jjj, iii] = CURRENT.STATE$log.p.theta + CURRENT.STATE$log.p.y.giv.theta.pseudo
##         mat.pseudo[jjj, iii] = CURRENT.STATE$log.p.theta
      }
  }
cat("\n")

APPROACH = "AA"
source("../functions/GP.FUNCTIONS.r")

mat.theta.giv.f = rep(0, length(psi.tau.vect))
CURRENT.STATE$psi.tau = TRUE.PSI.TAU
CURRENT.STATE$psi.sigma = log(TRUE.SIGMA)
CURRENT.STATE$f = TRUE.LATENT.F
for(iii in 1:length(psi.tau.vect))
  {
    CURRENT.STATE$psi.tau = psi.tau.vect[iii]

    compute.L.chol.Q.mat(CURRENT.STATE)
    compute.update.nu(CURRENT.STATE)
    compute.log.p.f.giv.theta(CURRENT.STATE)
    compute.log.p.theta(CURRENT.STATE)
    
    mat.theta.giv.f[iii] = CURRENT.STATE$log.p.theta + CURRENT.STATE$log.p.f.giv.theta
  }


mat.theta.giv.y.nu = rep(0, length(psi.tau.vect))
CURRENT.STATE$psi.tau = TRUE.PSI.TAU
CURRENT.STATE$psi.sigma = log(TRUE.SIGMA)
CURRENT.STATE$f = TRUE.LATENT.F
compute.L.chol.Q.mat(CURRENT.STATE)
compute.update.nu(CURRENT.STATE)
for(iii in 1:length(psi.tau.vect))
  {
    CURRENT.STATE$psi.tau = psi.tau.vect[iii]

    compute.utilities.t.giv.y.n.t(CURRENT.STATE)
    
    mat.theta.giv.y.nu[iii] = CURRENT.STATE$log.p.theta + CURRENT.STATE$log.p.y.giv.nu.theta
  }


APPROACH = "SURR"
source("../functions/GP.FUNCTIONS.r")

mat.theta.giv.g = matrix(0, length(psi.tau.vect), 1)

set.seed(123423)

for(kkk in 1:dim(mat.theta.giv.g)[2])
  {
    CURRENT.STATE$psi.tau = TRUE.PSI.TAU
    CURRENT.STATE$psi.sigma = log(TRUE.SIGMA)
    CURRENT.STATE$f = TRUE.LATENT.F

    compute.L.chol.Q.mat(CURRENT.STATE)
    compute.update.S.theta.R.theta(CURRENT.STATE)
    
    SURROGATE.G <<- sqrt(CURRENT.STATE$Stheta) * rnorm(DATA$n) + CURRENT.STATE$f

    compute.utilities.t.surrogate(CURRENT.STATE)

    for(iii in 1:length(psi.tau.vect))
      {
        CURRENT.STATE$psi.tau = psi.tau.vect[iii]

        compute.utilities.t.surrogate(CURRENT.STATE)
        
        mat.theta.giv.g[iii,kkk] = CURRENT.STATE$log.p.y.giv.f + CURRENT.STATE$log.p.g.giv.theta + CURRENT.STATE$log.p.theta
      }
}


mat.theta.giv.g = apply(mat.theta.giv.g, 1, mean)


lq.mat.pseudo = apply(mat.pseudo, 2, quantile, 0.025)
hq.mat.pseudo = apply(mat.pseudo, 2, quantile, 0.975)
median.mat.pseudo = apply(mat.pseudo, 2, quantile, 0.5)
mean.mat.pseudo = apply(mat.pseudo, 2, mean)

exp.mat.pseudo = exp(mat.pseudo - LOGSUM.VECT(mean.mat.pseudo) - log(step.psi.tau.vect))
exp.lq.mat.pseudo = apply(exp.mat.pseudo, 2, quantile, 0.025)
exp.hq.mat.pseudo = apply(exp.mat.pseudo, 2, quantile, 0.975)
exp.median.mat.pseudo = apply(exp.mat.pseudo, 2, quantile, 0.5)
exp.mean.mat.pseudo = apply(exp.mat.pseudo, 2, mean)

exp.mat.theta.giv.f = exp(mat.theta.giv.f - LOGSUM.VECT(mat.theta.giv.f) - log(step.psi.tau.vect))
exp.mat.theta.giv.y.nu = exp(mat.theta.giv.y.nu - LOGSUM.VECT(mat.theta.giv.y.nu) - log(step.psi.tau.vect))
exp.mat.theta.giv.g = exp(mat.theta.giv.g - LOGSUM.VECT(mat.theta.giv.g) - log(step.psi.tau.vect))


fileplot = paste("PLOT_COMPARE_GIBBS_n_", N, ".eps", sep="")

if(APPROXIMATION == "LA") name.APPROXIMATION = "LA"
if(APPROXIMATION == "EP") name.APPROXIMATION = "EP"

postscript(fileplot)
par("mar"=c(2.9,2.9,1.1,0.4), "las"=0, "mgp"=c(1.8,0.6,0))
## plot(exp(psi.tau.vect), exp.median.mat.pseudo, type="l", lwd=2, ylim = range(c(exp.lq.mat.pseudo,exp.hq.mat.pseudo,exp.mat.theta.giv.f)), ylab="posterior density", xlab="length-scale", col=1, main = paste("N = ", N, sep=""))
plot(exp(psi.tau.vect), exp.median.mat.pseudo, type="l", lwd=2, ylim = c(0,10), ylab="posterior density", xlab="length-scale", col=1, main = paste("N = ", N, sep=""))

points(exp(psi.tau.vect), exp.mat.theta.giv.f, type="l", lwd=2, lty=1, col=4)
points(exp(psi.tau.vect), exp.mat.theta.giv.y.nu, type="l", lwd=2, lty=1, col=2)
points(exp(psi.tau.vect), exp.mat.theta.giv.g, type="l", lwd=2, lty=1, col=3)
points(exp(psi.tau.vect), exp.median.mat.pseudo, type="l", lwd=2, lty=1, col=1)

abline(v=exp(TRUE.PSI.TAU), lty=2)

## legend(1.5, 25, c("SA", "AA", "SURR", "PM"), lwd=2, col=c(4, 2, 3, 1))
legend(1.5, 9, c("SA", "AA", "SURR", "PM"), lwd=2, col=c(4, 2, 3, 1))

dev.off()

