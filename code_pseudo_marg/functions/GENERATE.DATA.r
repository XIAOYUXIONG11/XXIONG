## Function to generate data from a GP probit model

GENERATE.DATA.GP = function(sigma, true.psi.tau, n, d, SEED.DATA=222)
  {
    x.gr = matrix(runif(2000*d), 2000, d)
    n.gr = dim(x.gr)[1]

    K.true = sigma * q.fun.xx(x.gr, sigma, exp(true.psi.tau))
    
    U.K.true = chol(K.true)
    set.seed(SEED.DATA)
    true.nu.vect = rnorm(n.gr)
    true.latent = t(U.K.true) %*% true.nu.vect

    p.y.gr = pnorm(true.latent)
    y.gr = rbinom(n.gr, 1, p.y.gr)*2 - 1

    ind.to.take = c(sample(which(y.gr == 1), n/2), sample(which(y.gr == -1), n/2))
    ind.shuffle = sample(c(1:n), n, replace=F)
    ind.to.take = ind.to.take[ind.shuffle]

    X = matrix(x.gr[ind.to.take,], ncol=d)
    y = y.gr[ind.to.take]

    K.true = sigma * q.fun.xx(X, sigma, exp(true.psi.tau))
    U.K.true = chol(K.true)
    zzz.tmp = solve(t(U.K.true)) %*% true.latent[ind.to.take]
    
    assign("TRUE.LATENT.NU", zzz.tmp, .GlobalEnv)
    assign("TRUE.LATENT.F", true.latent[ind.to.take], .GlobalEnv)
    assign("TEST.TRUE.LATENT", true.latent, .GlobalEnv)

    DATA$n = n
    DATA$d = d
    DATA$X = X
    DATA$y = y

    DATA$ntest = length(y.gr)
    DATA$X.TEST = x.gr
    DATA$y.TEST = y.gr
  }


## **************************************************
## ******************** DATASET
## **************************************************

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = D
set.seed(1)
TRUE.SIGMA = rgamma(1, shape=PRIOR.SIGMA.SHAPE, rate=PRIOR.SIGMA.RATE)
TRUE.PSI.TAU = runif(100, min=-1.5, max=-1)[1:NLS]
DATA = new.env()
GENERATE.DATA.GP(TRUE.SIGMA, TRUE.PSI.TAU, N, D, SEED.DATA=99090)
