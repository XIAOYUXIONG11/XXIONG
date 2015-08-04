## Function definitions for probit GP classification

reset.number.calls = function()
  {
    assign("NUMBER.CALLS.N3", 0, .GlobalEnv)
  }

logprior = function(psi.phi, psi.tau) 
  {
    dnorm(psi.phi, mean=PRIOR.PSI.PHI.MEAN, sd=sqrt(PRIOR.PSI.PHI.VAR), log=T) + PRIOR.TAU.SHAPE * psi.tau - PRIOR.TAU.RATE * exp(psi.tau)
  }

logistic = function(a)
  {
    1 / (1 + exp(-a))
  }

LOGSUM = function(a, b)
{
  if(min(a) >= min(b)) return(a + log(1 + exp(b-a)))
  b + log(1 + exp(a-b))
}

LOGDIFF = function(ab)
{
 ab[1] + log(1 - exp(ab[2]-ab[1]))
}

LOGABSDIFF = function(ab)
{
 if(ab[1] > ab[2]) return(ab[1] + log(1 - exp(ab[2]-ab[1])))
 ab[2] + log(1 - exp(ab[1]-ab[2]))
}

LOGSUM2 = function(ab)
{
  if(min(ab[1]) >= min(ab[2])) return(ab[1] + log(1 + exp(ab[2]-ab[1])))
  ab[2] + log(1 + exp(ab[1]-ab[2]))
}

LOGSUM.VECT = function(a)
{
  max.a = max(a)
  max.a + log(sum(exp(a - max.a)))
}

log.logistic = function(a)
  {
    -LOGSUM(-a, 0)
  }

## **************************************************************************************************** 
## **************************************************************************************************** 
## ************************************************** KERNEL FUNCTION
## Kernel function computing the nxn kernel matrix of n training points. The covariance is considered hyperspherical when w is a scalar or diagonal when w is a vector
## q.fun.xx = function(x, psi.phi, psi.tau)
##   {
##     phi = 2 * exp(log.logistic(psi.phi)) - 1
##     tau = exp(psi.tau)
    
##     res = matrix(0, nn+1, nn+1)

##     diag(res) = tau * (1 + phi^2)
##     res[1,1] = res[nn,nn] = tau
##     res[nn+1, nn+1] = PRIOR.MU.VAR

##     for(i in 1:nn) res[i, i-1] = res[i-1, i] = -phi * tau

##     as.spam(res)
##   }    


q.fun.xx = function(x, psi.phi, psi.tau)
  {
    nn = dim(x)[1]
    phi = 2 * exp(log.logistic(psi.phi)) - 1
    tau = exp(psi.tau)

    p1 = c(tau, -phi * tau)
    p2 = c(-phi * tau, tau * (1 + phi^2), -phi * tau)
    p3 = c(-phi * tau, tau, PRIOR.MU.VAR)
    zz = spam(0, nn+1, nn+1)
    zz@entries = c(p1, rep(p2, nn-2), p3)
    zz@colindices = as.integer(c(1,2, outer(c(1:3), c(0:(nn-2)), "+")))
    zz@rowpointers = as.integer(c(1, c(1:(nn-1)*3), (nn*3)-1, nn*3))

    zz
  }

compute.L.chol.Q.mat = function(STATE)
  {
    STATE$Q.mat.sparse = q.fun.xx(DATA$X, STATE$psi.phi, STATE$psi.tau)
    STATE$half.log.det.Q = sum(log(diag(chol(STATE$Q.mat.sparse))))
    STATE$sparsity.level = length(STATE$Q.mat.sparse@entries) / DATA$n^2
  }

## compute.log.p.y.giv.f = function(STATE)
##   {
##     STATE$log.p.y.giv.f = sum(dnorm(DATA$y, mean = 0, sd = exp((STATE$f + STATE$mu) / 2), log=T))
##     if(is.nan(STATE$log.p.y.giv.f) | is.infinite(STATE$log.p.y.giv.f)) STATE$log.p.y.giv.f = -Inf
##   }

## compute.logjointlik.AR.f.giv.y.f.t = function(STATE)
## {
##   STATE$log.p.y.giv.f + STATE$log.p.f.giv.theta + STATE$log.p.theta
## }

compute.log.p.theta = function(STATE)
{
  STATE$log.p.theta = logprior(STATE$psi.phi, STATE$psi.tau)
}


## **************************************************************************************************** PM

compute.initial.utilities = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.log.p.y.giv.theta.pseudo.marginal(STATE)
  compute.log.p.theta(STATE)
  
  compute.logjointlik(STATE)
}

compute.utilities.pseudo.marginal = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.log.p.y.giv.theta.pseudo.marginal(STATE)
  compute.log.p.theta(STATE)
}

fun.log.p.y.giv.h = function(h)
  {
    f = h[1:DATA$n] + h[DATA$n+1]
    res = sum(dnorm(DATA$y, mean = 0, sd = exp((f) / 2), log=T))
    
    if(is.nan(res)) res = -Inf
    res
  }

fun.log.p.h.giv.theta = function(h, STATE)
{
  Q.h = STATE$Q.mat.sparse %*% h
  STATE$half.log.det.Q - 0.5 * crossprod(h, Q.h)
}

compute.log.p.y.giv.theta.pseudo.marginal = function(STATE)
{
  if(APPROXIMATION == "LA") approx.model = LA.sparse(STATE)

  h.pseudo = t(rmvnorm.prec(NPSEUDO, approx.model$mu, approx.model$P))

  tmp = approx.model$P %*% (h.pseudo - approx.model$mu)
  crossprodtmp = apply((h.pseudo - approx.model$mu) * tmp, 2, sum)

  log.det.P = sum(log(diag(chol(approx.model$P))))
  
  log.pseudo.weights = apply(h.pseudo, 2, fun.log.p.y.giv.h) + apply(h.pseudo, 2, fun.log.p.h.giv.theta, STATE) - log.det.P + 0.5 * crossprodtmp - log(NPSEUDO)
  
  res = LOGSUM.VECT(log.pseudo.weights)
  if(is.nan(res)) res = -Inf
  
  STATE$log.p.y.giv.theta.pseudo = res
}

compute.logjointlik = function(STATE)
{
  STATE$logjointlik = STATE$log.p.y.giv.theta.pseudo + STATE$log.p.theta

  for(i in 1:SAMPLER$nblocks)
    STATE[[paste("LOGJOINTLIK.AR", i, sep="")]] = SAMPLER$blocks[[i]]$compute.logjointlik.AR(STATE)
}



compute.logjointlik.AR.pseudo.marginal = function(STATE)
{
  STATE$log.p.y.giv.theta.pseudo + STATE$log.p.theta
}
