## Code to generate the results to produce figure 2
## Assessment of the importance sampling distribution

## **************************************************
## ******************** Initial parameters
## **************************************************

## Remove comments according to N, APPROXIMATION and NPSEUDO (N_imp in the paper)

## N = 10
## NPSEUDO = 1

## APPROXIMATION = "LA"
## APPROXIMATION = "EP"


COVARIANCE.MODE = "H"
NREP = 500
D = 2

## **************************************************
## ******************** Load useful functions
## **************************************************

APPROACH = "PM"

pdf.options(width=10,height=10,pointsize=32)
ps.options(width=10,height=10,paper="special",horizontal=F, pointsize=32)
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
max.tau.vect = 2.8

psi.tau.vect = log(seq(min.tau.vect, max.tau.vect, length.out=40))
step.psi.tau.vect = diff(seq(min.tau.vect, max.tau.vect, length.out=40)[c(1,2)])
mat.pseudo = matrix(0, NREP, length(psi.tau.vect))

for(iii in 1:length(psi.tau.vect))
  {
    cat(formatC(iii/length(psi.tau.vect)*100, format="f", digits=2), "%    \r")
    CURRENT.STATE$psi.tau = psi.tau.vect[iii]
    CURRENT.STATE$psi.sigma = log(TRUE.SIGMA)

    compute.L.chol.Q.mat(CURRENT.STATE)
    compute.log.p.theta(CURRENT.STATE)

    for(jjj in 1:NREP)
      {
        compute.log.p.y.giv.theta.pseudo.marginal(CURRENT.STATE)

        mat.pseudo[jjj, iii] = CURRENT.STATE$log.p.theta + CURRENT.STATE$log.p.y.giv.theta.pseudo
      }
  }
cat("\n")

lq.mat.pseudo = apply(mat.pseudo, 2, quantile, 0.025)
hq.mat.pseudo = apply(mat.pseudo, 2, quantile, 0.975)
median.mat.pseudo = apply(mat.pseudo, 2, quantile, 0.5)
mean.mat.pseudo = apply(mat.pseudo, 2, mean)

exp.mat.pseudo = exp(mat.pseudo - LOGSUM.VECT(mean.mat.pseudo) - log(step.psi.tau.vect))
exp.lq.mat.pseudo = apply(exp.mat.pseudo, 2, quantile, 0.025)
exp.hq.mat.pseudo = apply(exp.mat.pseudo, 2, quantile, 0.975)
exp.median.mat.pseudo = apply(exp.mat.pseudo, 2, quantile, 0.5)
exp.mean.mat.pseudo = apply(exp.mat.pseudo, 2, mean)

filesave = paste("RES_PLOT_", "PROBIT", "_n_", N, "_method_", APPROXIMATION, "_nimp_", NPSEUDO, ".Rdata", sep="")

save.image(filesave)

