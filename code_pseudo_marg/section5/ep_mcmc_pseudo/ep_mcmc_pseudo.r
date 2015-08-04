## PM approach to sample hyper-parameters and Elliptical Slice Sampling for latent variables

source("../../functions/batchmeans.R")

APPROACH = "PM"
source("../../functions/GENERAL.SAMPLER.REAL.DATA.r")
source("../../functions/GP.FUNCTIONS.r")
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


OMEGA = 1e-6

PRIOR.PAR.SHAPE = 1
PRIOR.PAR.RATE = 0.1 ## / sqrt(D)

PRIOR.SIGMA.SHAPE = 1.2
PRIOR.SIGMA.RATE = 0.2

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = D

SCHEME.SAMPLER.HYPER = "MH"
SEED = 101
WITH.PREDICTIONS = F

APPROXIMATION = "EP"

NPSEUDO = 64

CHECK.MODE = F

BLOCKING = list()

BLOCKING[[1]] = list(group=c("psi.sigma", "psi.tau"), scheme="pseudo.marginal", method=SCHEME.SAMPLER.HYPER)
BLOCKING[[2]] = list(group="f", scheme="f.giv.y.f.t", method="ELLIPTICALSS")

reset.number.calls()


## **************************************************
## ******************** MODEL
## **************************************************

compute.approximation = EP.NOGRAD

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
SAMPLER$nsamples = 5000
SAMPLER$nburnin = 1000
SAMPLER$nbatch = 100
SAMPLER$save.every = 1

set.seed(SEED)
SAMPLER$init.par$psi.tau = log(rgamma(MODEL$NLS, shape=PRIOR.PAR.SHAPE, rate=PRIOR.PAR.RATE))
SAMPLER$init.par$psi.sigma = log(rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE))
SAMPLER$init.par$f = rep(0, DATA$n)

ADAPTIVE = F

MAXTOSTORE = min(SAMPLER$nsamples, 5000)

FILESAVECHAINS.psi.tau = paste("CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.psi.sigma = paste("CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
FILESAVECHAINS.f = paste("CHAINS_F_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")

set.seed(SEED)
res = GENERAL.SAMPLER()

TEST.PREDICTIONS.MCMC = matrix(0, DATA$ntest, 0)

for(iii in (SAMPLER$nburnin+1):SAMPLER$nsamples)
  {
    psi.sigma = SAMPLES$psi.sigma[iii]
    psi.tau = SAMPLES$psi.tau[iii]
    f = SAMPLES$f[iii,]
    
    L.chol = sqrt(exp(psi.sigma)) * t(chol(q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))))

    K.star.star = rep(0, DATA$ntest)
    K.star = exp(psi.sigma) * q.fun.xy(DATA$X, DATA$X.TEST, exp(psi.tau))
    for(i in 1:DATA$ntest) K.star.star[i] = exp(psi.sigma) * q.fun.xx(matrix(DATA$X.TEST[i,], nrow=1), exp(psi.sigma), exp(psi.tau))

    Lf = forwardsolve(L.chol, f)
    LKstar = forwardsolve(L.chol, K.star)
    
    mstar = t(LKstar) %*% Lf ## t(K.star) %*% solve(K) %*% f
    vstar = K.star.star - apply(LKstar^2, 2, sum) ## K.star.star - apply((solve(K) %*% K.star) * K.star, 2, sum)
    vstar[vstar < 0] = 1e-9
    
    fstar = matrix(rnorm(DATA$ntest*100, mean=mstar, sd=sqrt(vstar)), nrow=DATA$ntest)

    p.y.giv.f.star = apply(exp(pnorm(fstar, log=T)), 1, mean)
    
    TEST.PREDICTIONS.MCMC = cbind(TEST.PREDICTIONS.MCMC, p.y.giv.f.star)
  }

tmp = bmmat(t(TEST.PREDICTIONS.MCMC))
row.names(tmp) = NULL

predictions.test = tmp[,1]
predictions.test.se = tmp[,2]
