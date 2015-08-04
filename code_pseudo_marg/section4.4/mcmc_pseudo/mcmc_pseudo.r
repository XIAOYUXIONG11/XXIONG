## PM approach to sample hyper-parameters 

APPROACH = "PM"
NPSEUDO = 1
## SEED = 1

source("../../functions/batchmeans.R")

source("../../functions/GENERAL.SAMPLER.SECTION.4.4.r")
source("../../functions/GP.FUNCTIONS.r")
source("../../functions/LA.r")
source("../../functions/EP.r")
source("../../functions/MH.r")
source("../../functions/ELLIPTICALSS.r")


compute.utilities.f.giv.y.f.t = function(STATE)
{
  compute.update.nu(STATE)
  compute.log.p.y.giv.f(STATE)
  compute.log.p.y.giv.nu.theta(STATE)
  compute.log.p.nu(STATE)
  compute.log.p.f.giv.theta(STATE)
}


compute.update.nu = function(STATE)
{
  STATE$nu = forwardsolve(STATE$L.chol.Q.mat, STATE$f) / exp(STATE$psi.sigma/2)
}

compute.log.p.y.giv.nu.theta = function(STATE)
  {
    STATE$log.p.y.giv.nu.theta = sum(pnorm(DATA$y * STATE$f, log=T)) ## sum(DATA$y * STATE$log.logistic.vect + (1 - DATA$y) * STATE$log.logistic.vect.minus)
    if(is.nan(STATE$log.p.y.giv.nu.theta) | is.infinite(STATE$log.p.y.giv.nu.theta)) STATE$log.p.y.giv.nu.theta = -Inf
  }

compute.log.p.nu = function(STATE)
{
  STATE$log.p.nu = - 0.5 * crossprod(STATE$nu)
}

compute.log.p.f.giv.theta = function(STATE)
{
  STATE$log.p.f.giv.theta = - sum(log(diag(STATE$L.chol.Q.mat))) - DATA$n / 2 * STATE$psi.sigma - 0.5 * crossprod(STATE$nu)
}

compute.initial.utilities = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.utilities.f.giv.y.f.t(STATE)
  
  compute.log.p.y.giv.theta.pseudo.marginal(STATE)
  compute.log.p.theta(STATE)
  
  compute.logjointlik(STATE)
}


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

BLOCKING[[1]] = list(group=c("psi.sigma", "psi.tau"), scheme="pseudo.marginal", method=SCHEME.SAMPLER.HYPER, repetitions=1)

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

SAMPLER$blocks[[1]]$cov.proposal = diag(c(0.5, rep(0.5, MODEL$NLS)))
SAMPLER$blocks[[1]]$L.chol.cov.proposal = t(chol(SAMPLER$blocks[[1]]$cov.proposal))

## ***************************************************
## ******************** SAMPLING FROM A RANDOM PLACE DEPENDING ON THE SEED 
## ***************************************************
compute.log.p.y.giv.theta.pseudo.marginal = function(STATE)
{
##   approx.model <<- EP.NOGRAD.CONVERGENCE.ON.Z(STATE$psi.sigma, STATE$psi.tau)
  approx.model <<- LA.Z(STATE)
  
  STATE$log.p.y.giv.theta.pseudo = approx.model$log.Z.LA
}
  
##   tmp = forwardsolve(approx.model$L.Sigma, f.pseudo - c(approx.model$mu))
##   crossprodtmp = apply(tmp^2, 2, sum)
##   log.pseudo.weights = apply(f.pseudo, 2, fun.log.p.y.giv.f) + apply(f.pseudo, 2, fun.log.p.f.giv.theta, STATE) + sum(log(diag(approx.model$L.Sigma))) + 0.5 * crossprodtmp - log(NPSEUDO)


SAMPLER$nsamples = 2000
SAMPLER$nburnin = 2000
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

set.seed(SEED+876)
SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))
SAMPLER$init.par$f = rep(0, DATA$n)

ADAPTIVE = T

MAXTOSTORE = min(SAMPLER$nsamples, 2000)

FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.f = paste("CHAINS_F_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")

set.seed(SEED+876)
time.LA = system.time(assign("res", GENERAL.SAMPLER(), envir=.GlobalEnv))[3]

## ***************************************************
## ******************** SAMPLING FROM A RANDOM PLACE DEPENDING ON THE SEED 
## ***************************************************
source("../../functions/GP.FUNCTIONS.r")

SAMPLER$init.par$psi.tau = SAMPLES$psi.tau[SAMPLER$nsamples,]
SAMPLER$init.par$psi.sigma = SAMPLES$psi.sigma[SAMPLER$nsamples]

SAMPLER$init.par$f = SAMPLES$f[SAMPLER$nsamples,]

SAMPLER$nsamples = 10000
SAMPLER$nburnin = 0
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

## set.seed(SEED)
## ADAPTIVE = T

MAXTOSTORE = min(SAMPLER$nsamples, 2000)
## set.seed(SEED)
time.PM = system.time(assign("res", GENERAL.SAMPLER(), envir=.GlobalEnv))[3]

FILESAVE.times = paste("TIMES_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
cat(time.LA, "\n", time.PM, "\n", file=FILESAVE.times)
