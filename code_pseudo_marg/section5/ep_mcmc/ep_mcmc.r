## MCMC inference of the hyper-parameters in probit GP classification.
## The sampling is based on the marginal likelihood obtained using EP

source("../../functions/batchmeans.R")

APPROACH = "PM"
source("../../functions/GENERAL.SAMPLER.REAL.DATA.r")
source("../../functions/GP.FUNCTIONS.r")
source("../../functions/EP.r")
source("../../functions/MH.r")

compute.log.p.y.giv.theta.pseudo.marginal = function(STATE)
  {
    tmp <<- EP.NOGRAD(STATE$psi.sigma, STATE$psi.tau)
    STATE$log.p.y.giv.theta.pseudo = tmp$log.Z.ep

    if(WITH.PREDICTIONS == T) STATE$EP.MODEL = tmp
  }

OMEGA = 1e-6

PRIOR.PAR.SHAPE = 1
PRIOR.PAR.RATE = 0.1 ## / sqrt(D)

PRIOR.SIGMA.SHAPE = 1.2
PRIOR.SIGMA.RATE = 0.2

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = D


SCHEME.SAMPLER.HYPER = "MH"
SEED = 101
WITH.PREDICTIONS = T

CHECK.MODE = F

BLOCKING = list()

BLOCKING[[1]] = list(group=c("psi.sigma", "psi.tau"), scheme="pseudo.marginal", method=SCHEME.SAMPLER.HYPER)

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

SAMPLER$blocks[[1]]$cov.proposal = diag(c(0.5, rep(0.5, MODEL$NLS)))
SAMPLER$blocks[[1]]$L.chol.cov.proposal = t(chol(SAMPLER$blocks[[1]]$cov.proposal))

## ***************************************************
## ******************** SAMPLING FROM A RANDOM PLACE DEPENDING ON THE SEED (ADAPT FOR A SMALL NUMBER OF STEPS JUST TO MOVE OUT FROM THE STARTING POSITION IF THE CHAIN GETS STUCK)
## ***************************************************
SAMPLER$nsamples = 5000
SAMPLER$nburnin = 1000
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

set.seed(SEED)
SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))

ADAPTIVE = F

MAXTOSTORE = min(SAMPLER$nsamples, 5000)

FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")

set.seed(SEED)
res = GENERAL.SAMPLER()

tmp = bmmat(t(TEST.PREDICTIONS.MCMC))
row.names(tmp) = NULL

predictions.test = tmp[,1]
predictions.test.se = tmp[,2]
