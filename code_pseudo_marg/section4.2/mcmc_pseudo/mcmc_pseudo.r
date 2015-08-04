## Pseudo marginal approach

source("../../functions/EP.r")
source("../../functions/LA.r")
source("../../functions/GENERAL.SAMPLER.r")
source("../../functions/MH.r")
source("../../functions/ELLIPTICALSS.r")

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = DTRAINING

SCHEME.SAMPLER.HYPER = "MH"
WITH.PREDICTIONS = F

if(APPROXIMATION == "EP") PSEUDO.TYPE = "EP"
if(APPROXIMATION == "LA") PSEUDO.TYPE = "LA"

CHECK.MODE = F

BLOCKING = list()

BLOCKING[[1]] = list(group=c("psi.sigma", "psi.tau"), scheme="pseudo.marginal", method=SCHEME.SAMPLER.HYPER, repetitions=1)

reset.number.calls()

## **************************************************
## ******************** MODEL
## **************************************************

MODEL = new.env()
MODEL$NLS = NLS
MODEL$SIZE.GROUPS = c(1, MODEL$NLS)
names(MODEL$SIZE.GROUPS) = MODEL$NAMES.GROUPS = c("psi.sigma", "psi.tau")

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

## ******************** MH for jointly update psi.sigma and psi.tau
if((NTRAINING == 50) & (DTRAINING==2) & (COVARIANCE.MODE == "H")) SAMPLER$blocks[[1]]$cov.proposal = diag(2, MODEL$NLS + 1)
if((NTRAINING == 50) & (DTRAINING==10) & (COVARIANCE.MODE == "H")) SAMPLER$blocks[[1]]$cov.proposal = diag(5, MODEL$NLS + 1)
if((NTRAINING == 200) & (DTRAINING==2) & (COVARIANCE.MODE == "H")) SAMPLER$blocks[[1]]$cov.proposal = diag(1, MODEL$NLS + 1)
if((NTRAINING == 200) & (DTRAINING==10) & (COVARIANCE.MODE == "H")) SAMPLER$blocks[[1]]$cov.proposal = diag(1, MODEL$NLS + 1)

if((NTRAINING == 50) & (DTRAINING==2) & (COVARIANCE.MODE == "D")) SAMPLER$blocks[[1]]$cov.proposal = diag(2, MODEL$NLS + 1)
if((NTRAINING == 50) & (DTRAINING==10) & (COVARIANCE.MODE == "D")) SAMPLER$blocks[[1]]$cov.proposal = diag(0.5, MODEL$NLS + 1)
if((NTRAINING == 200) & (DTRAINING==2) & (COVARIANCE.MODE == "D")) SAMPLER$blocks[[1]]$cov.proposal = diag(0.2, MODEL$NLS + 1)
if((NTRAINING == 200) & (DTRAINING==10) & (COVARIANCE.MODE == "D")) SAMPLER$blocks[[1]]$cov.proposal = diag(0.5, MODEL$NLS + 1)

SAMPLER$blocks[[1]]$L.chol.cov.proposal = t(chol(SAMPLER$blocks[[1]]$cov.proposal))

## ***************************************************
## ******************** SAMPLING FROM A RANDOM PLACE SAMPLED FROM THE PRIOR DEPENDING ON THE SEED
## ***************************************************
SAMPLER$nsamples = 15000
SAMPLER$nburnin = 5000
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

set.seed(SEED)
SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))

MAXTOSTORE = min(SAMPLER$nsamples, 15000)

FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_", NTRAINING, "_", DTRAINING, "_", APPROXIMATION, "_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_", NTRAINING, "_", DTRAINING, "_", APPROXIMATION, "_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVE.NON3 = paste("NON3_", NTRAINING, "_", DTRAINING, "_", APPROXIMATION, "_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")

set.seed(SEED)
res = GENERAL.SAMPLER()
