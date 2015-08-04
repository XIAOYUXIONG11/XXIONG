## Surrogate method

source("../../functions/batchmeans.R")

APPROACH = "SURR"
source("../../functions/GENERAL.SAMPLER.SECTION.4.4.r")
source("../../functions/GP.FUNCTIONS.SECTION.4.4.r")
source("../../functions/LA.r")
source("../../functions/MH.r")
source("../../functions/ELLIPTICALSS.r")

OMEGA = 1e-9
if(DATASET == "breast_cancer_wisconsin") OMEGA = 1e-8

PRIOR.PAR.SHAPE = 1
if(COVARIANCE.MODE == "H") PRIOR.PAR.RATE = 1 / sqrt(DATA$d)
if(COVARIANCE.MODE == "D") PRIOR.PAR.RATE = 1

PRIOR.SIGMA.SHAPE = 1.1
PRIOR.SIGMA.RATE = 0.1

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = DATA$d

SCHEME.SAMPLER.HYPER = "MH"
WITH.PREDICTIONS = F

APPROXIMATION = "LA"

CHECK.MODE = F

BLOCKING = list()

BLOCKING[[1]] = list(group="f", scheme="f.giv.y.f.t", method="ELLIPTICALSS", repetitions=10)
BLOCKING[[2]] = list(group=c("psi.sigma", "psi.tau"), scheme="t.surrogate", method=SCHEME.SAMPLER.HYPER, repetitions=1)

reset.number.calls()


## **************************************************
## ******************** MODEL
## **************************************************

compute.approximation = LA

MODEL = new.env()
MODEL$NLS = NLS
MODEL$SIZE.GROUPS = c(1, MODEL$NLS, DATA$n)
names(MODEL$SIZE.GROUPS) = MODEL$NAMES.GROUPS = c("psi.sigma", "psi.tau", "f")

FILL.SAMPLER = function(SAMPLER)
{
  SAMPLER$blocks = list()
  SAMPLER$nblocks = length(BLOCKING)

  for(ind.block in 1:SAMPLER$nblocks)
    {
      SAMPLER$blocks[[ind.block]] = list()
      SAMPLER$blocks[[ind.block]]$group = BLOCKING[[ind.block]]$group
      SAMPLER$blocks[[ind.block]]$scheme = BLOCKING[[ind.block]]$scheme
      SAMPLER$blocks[[ind.block]]$method = BLOCKING[[ind.block]]$method
      SAMPLER$blocks[[ind.block]]$repetitions = BLOCKING[[ind.block]]$repetitions

      SAMPLER$blocks[[ind.block]]$compute.utilities = get(paste("compute.utilities.", SAMPLER$blocks[[ind.block]]$scheme, sep=""))
      SAMPLER$blocks[[ind.block]]$compute.logjointlik.AR = get(paste("compute.logjointlik.AR.", SAMPLER$blocks[[ind.block]]$scheme, sep=""))

      SAMPLER$blocks[[ind.block]]$names.logjointlik.AR = paste("logjointlik.AR.", paste(BLOCKING[[ind.block]]$group, collapse="."), sep="")
    }
  
  SAMPLER
}

SAMPLER = new.env()
SAMPLER$init.par = list()

SAMPLER$compute.initial.utilities = compute.initial.utilities
SAMPLER$compute.logjointlik = compute.logjointlik
FILL.SAMPLER(SAMPLER)

SAMPLER$blocks[[2]]$cov.proposal = diag(0.05, MODEL$NLS + 1)
SAMPLER$blocks[[2]]$L.chol.cov.proposal = t(chol(SAMPLER$blocks[[2]]$cov.proposal))

## ***************************************************
## ******************** SAMPLING FROM A RANDOM PLACE DEPENDING ON THE SEED (ADAPT FOR A SMALL NUMBER OF STEPS JUST TO MOVE OUT FROM THE STARTING POSITION IF THE CHAIN GETS STUCK)
## ***************************************************
SAMPLER$nsamples = 12000
SAMPLER$nburnin = 2000
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

set.seed(SEED+876)
SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))
SAMPLER$init.par$f = rep(0, DATA$n)
## SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
## SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))
## SAMPLER$init.par$nu = rnorm(DATA$n)
## SAMPLER$init.par$f = sqrt(exp(SAMPLER$init.par$psi.sigma)) * t(chol(q.fun.xx(DATA$X, exp(SAMPLER$init.par$psi.sigma), exp(SAMPLER$init.par$psi.tau)))) %*% SAMPLER$init.par$nu

ADAPTIVE = T

MAXTOSTORE = min(SAMPLER$nsamples, 2000)

FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.f = paste("CHAINS_F_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
## FILESAVECHAINS.f = paste("CHAINS_F_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")
## FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")
## FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")
## FILESAVE.NON3 = paste("NON3_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")

## set.seed(SEED)
time.SURR = system.time(assign("res", GENERAL.SAMPLER(), envir=.GlobalEnv))[3]

FILESAVE.times = paste("TIMES_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
cat(time.SURR, "\n", file=FILESAVE.times)
