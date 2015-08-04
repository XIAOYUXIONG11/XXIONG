## Ancillary Aumentation sampling scheme
## 

## **************************************************
## ******************** Initial parameters
## **************************************************

SCHEME.SAMPLER.LATENT = "ELLIPTICALSS"
SCHEME.SAMPLER.HYPER = "MH"

## **************************************************
## ******************** Load useful functions
## **************************************************

## ps.options(width=12,height=9,paper="special",horizontal=F, pointsize=24)
pdf.options(width=10,height=10,pointsize=24)

source("../../functions/GENERAL.SAMPLER.r")
source("../../functions/MH.r")
source("../../functions/ELLIPTICALSS.r")

CHECK.MODE = F

BLOCKING = list()
BLOCKING[[1]] = list(group="f", scheme="f.giv.y.f.t", method=SCHEME.SAMPLER.LATENT, repetitions=10)
BLOCKING[[2]] = list(group=c("psi.sigma", "psi.tau"), scheme="t.giv.y.n.t", method=SCHEME.SAMPLER.HYPER, repetitions=1)

reset.number.calls()

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = DTRAINING

MODEL = new.env()
MODEL$NLS = NLS
MODEL$SIZE.GROUPS = c(DATA$n, DATA$n, 1, MODEL$NLS)
names(MODEL$SIZE.GROUPS) = MODEL$NAMES.GROUPS = c("f", "nu", "psi.sigma", "psi.tau")

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
    }
  
  SAMPLER
}

SAMPLER = new.env()
SAMPLER$init.par = list()

SAMPLER$compute.initial.utilities = compute.initial.utilities
SAMPLER$compute.logjointlik = compute.logjointlik
FILL.SAMPLER(SAMPLER)


## [2] ******************** MH with AA for jointly update psi.sigma and psi.tau
SAMPLER$blocks[[2]]$cov.proposal = diag(0.2, MODEL$NLS + 1)
SAMPLER$blocks[[2]]$L.chol.cov.proposal = t(chol(SAMPLER$blocks[[2]]$cov.proposal))


## ***************************************************
## ******************** SAMPLING FROM A RANDOM PLACE DEPENDING ON THE SEED (ADAPT FOR A SMALL NUMBER OF STEPS JUST TO MOVE OUT FROM THE STARTING POSITION IF THE CHAIN GETS STUCK)
## ***************************************************
SAMPLER$nsamples = 15000
SAMPLER$nburnin = 5000
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

set.seed(SEED)
SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))
SAMPLER$init.par$nu = rnorm(DATA$n)
SAMPLER$init.par$f = sqrt(exp(SAMPLER$init.par$psi.sigma)) * t(chol(q.fun.xx(DATA$X, exp(SAMPLER$init.par$psi.sigma), exp(SAMPLER$init.par$psi.tau)))) %*% SAMPLER$init.par$nu

MAXTOSTORE = min(SAMPLER$nsamples, 15000)

FILESAVECHAINS.f = paste("CHAINS_F_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")
FILESAVE.NON3 = paste("NON3_COV_", COVARIANCE.MODE, "_N_", DATA$n, "_D_", DATA$d, "_SEED_", SEED, ".txt", sep="")

set.seed(SEED)
ELAPSED.TIME = system.time(assign("res", GENERAL.SAMPLER()))

rm(CURRENT.STATE, PROPOSED.STATE, SAMPLER)
rm(nu, envir=SAMPLES)

