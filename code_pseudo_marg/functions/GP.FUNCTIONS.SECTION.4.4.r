## Function definitions for probit GP classification

reset.number.calls = function()
  {
    assign("NUMBER.CALLS.N3", 0, .GlobalEnv)
  }

logprior = function(psi.sigma, psi.tau) 
  {
    PRIOR.SIGMA.SHAPE * psi.sigma - PRIOR.SIGMA.RATE * exp(psi.sigma) + sum(PRIOR.PAR.SHAPE * psi.tau - PRIOR.PAR.RATE * exp(psi.tau))
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
q.fun.xx = function(x, sigma, w)
  {
    d = dim(x)[2]
    n = dim(x)[1]
    w[w < 1e-19] = 1e-19
    W = matrix(0, d, d); diag(W) = 1 / w^2

    tmp0 = x %*% W
    tmp1 = apply(tmp0 * x, 1, sum)
    tmp3 = tmp0 %*% t(x)
    
    distances2 = outer(tmp1, tmp1, "+") - 2 * tmp3
    distances2[distances2 < 0] = 0
    
    exp(-0.5 * distances2) + diag(OMEGA/sigma, n)
  }

q.fun.xy = function(x, y, w)
  {
    d = dim(x)[2]
    w[w < 1e-19] = 1e-19
    W = matrix(0, d, d); diag(W) = 1 / w^2

    tmp1 = apply( (x %*% W) * x, 1, sum )
    tmp2 = apply( (y %*% W) * y, 1, sum )
    tmp3 = x %*% W %*% t(y)

    distances2 = outer(tmp1, tmp2, "+") - 2 * tmp3

    exp( -0.5 * distances2 )
  }

compute.L.chol.Q.mat = function(STATE)
  {
    STATE$Q.mat = q.fun.xx(DATA$X, exp(STATE$psi.sigma), exp(STATE$psi.tau))

    STATE$L.chol.Q.mat = t(chol(STATE$Q.mat))
    assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+1, envir=.GlobalEnv)
  }

compute.log.p.y.giv.f = function(STATE)
  {
    STATE$log.p.y.giv.f = sum(pnorm(DATA$y * STATE$f, log=T)) ## sum(DATA$y * STATE$log.logistic.vect + (1 - DATA$y) * STATE$log.logistic.vect.minus)
    if(is.nan(STATE$log.p.y.giv.f) | is.infinite(STATE$log.p.y.giv.f)) STATE$log.p.y.giv.f = -Inf
  }

compute.logjointlik.AR.f.giv.y.f.t = function(STATE)
{
  STATE$log.p.y.giv.f + STATE$log.p.f.giv.theta + STATE$log.p.theta
}

compute.log.p.theta = function(STATE)
{
  STATE$log.p.theta = logprior(STATE$psi.sigma, STATE$psi.tau)
}

## **************************************************************************************************** AA
if(APPROACH == "AA") {
compute.update.nu = function(STATE)
{
  STATE$nu = forwardsolve(STATE$L.chol.Q.mat, STATE$f) / exp(STATE$psi.sigma/2)
}

compute.update.f = function(STATE)
{
  STATE$f = STATE$L.chol.Q.mat %*% STATE$nu * exp(STATE$psi.sigma/2)
}

compute.log.p.y.giv.nu.theta = function(STATE)
  {
    STATE$log.p.y.giv.nu.theta = sum(pnorm(DATA$y * STATE$f, log=T)) ## sum(DATA$y * STATE$log.logistic.vect + (1 - DATA$y) * STATE$log.logistic.vect.minus)
    if(is.nan(STATE$log.p.y.giv.nu.theta) | is.infinite(STATE$log.p.y.giv.nu.theta)) STATE$log.p.y.giv.nu.theta = -Inf
  }

compute.log.p.f.giv.theta = function(STATE)
{
  STATE$log.p.f.giv.theta = - sum(log(diag(STATE$L.chol.Q.mat))) - DATA$n / 2 * STATE$psi.sigma - 0.5 * crossprod(STATE$nu)
}

compute.log.p.nu = function(STATE)
{
  STATE$log.p.nu = - 0.5 * crossprod(STATE$nu)
}

compute.initial.utilities = function(STATE)
{
  compute.L.chol.Q.mat(STATE)

  compute.update.nu(STATE)
  
  compute.log.p.y.giv.f(STATE)
  compute.log.p.y.giv.nu.theta(STATE)
  compute.log.p.f.giv.theta(STATE)
  compute.log.p.nu(STATE)
  compute.log.p.theta(STATE)
  compute.logjointlik(STATE)
}

compute.utilities.f.giv.y.f.t = function(STATE)
{
  compute.update.nu(STATE)
  compute.log.p.y.giv.f(STATE)
  compute.log.p.y.giv.nu.theta(STATE)
  compute.log.p.nu(STATE)
  compute.log.p.f.giv.theta(STATE)
}

compute.utilities.t.giv.y.n.t = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.update.f(STATE)

  compute.log.p.y.giv.nu.theta(STATE)
  compute.log.p.y.giv.f(STATE)
  compute.log.p.f.giv.theta(STATE)
  compute.log.p.theta(STATE)
}

compute.logjointlik = function(STATE)
{
  STATE$logjointlik = STATE$log.p.y.giv.f + STATE$log.p.f.giv.theta + STATE$log.p.theta

  for(i in 1:SAMPLER$nblocks)
    STATE[[paste("LOGJOINTLIK.AR", i, sep="")]] = SAMPLER$blocks[[i]]$compute.logjointlik.AR(STATE)
}

compute.logjointlik.AR.t.giv.y.n.t = function(STATE)
{
  STATE$log.p.y.giv.nu.theta + STATE$log.p.nu + STATE$log.p.theta
}


}

## **************************************************************************************************** SURR
if(APPROACH == "SURR") {

compute.update.mtg = function(STATE)
{
  STATE$mtg = STATE$Rtheta %*% (SURROGATE.G / STATE$Stheta)
}

compute.update.nu = function(STATE)
{
  STATE$nu = forwardsolve(STATE$L.chol.R.mat, STATE$f - STATE$mtg)
}

compute.update.f = function(STATE)
{
  STATE$f = STATE$L.chol.R.mat %*% STATE$nu + STATE$mtg
}

compute.log.p.f.giv.theta = function(STATE)
{
  zzz = forwardsolve(STATE$L.chol.Q.mat, STATE$f) / exp(STATE$psi.sigma/2)
  STATE$log.p.f.giv.theta = - sum(log(diag(STATE$L.chol.Q.mat))) - DATA$n / 2 * STATE$psi.sigma - 0.5 * crossprod(zzz)
}

compute.update.S.theta.R.theta = function(STATE)
{
  prior_var = exp(STATE$psi.sigma)
  prior_precision = 1/prior_var
  ## post_var = prior_var * (1 - 1/(pi/2 + 4/prior_var)) ## LOGISTIC REGRESSION
  post_var = prior_var * (1 - 1/(pi/2 + pi/2/prior_var)) ## PROBIT REGRESSION
  
  post_precision = 1/post_var;
  aux_var = rep(1 / (post_precision - prior_precision), DATA$n)
  
  STATE$Stheta = c(aux_var)
##   STATE$Rtheta = solve(exp(-STATE$psi.sigma) * solve(STATE$Q.mat) + diag(1/STATE$Stheta))
  STATE$Rtheta = diag(STATE$Stheta) - t(t((STATE$Stheta) * solve(diag(STATE$Stheta) + STATE$Q.mat * exp(STATE$psi.sigma))) * (STATE$Stheta))
  STATE$L.chol.R.mat = t(chol(STATE$Rtheta))

  assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+2, envir=.GlobalEnv) ## Inversion to compute R and Cholesky of R
}

compute.log.p.g.giv.theta = function(STATE)
{
  tmp = t(chol(STATE$Q.mat * exp(STATE$psi.sigma) + diag(STATE$Stheta)))
  STATE$log.p.g.giv.theta = - sum(log(diag(tmp))) - 0.5 * crossprod(forwardsolve(tmp, SURROGATE.G))
}

compute.initial.utilities = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.update.S.theta.R.theta(STATE)

  compute.log.p.y.giv.f(STATE)
  compute.log.p.f.giv.theta(STATE)
  compute.log.p.theta(STATE)
  compute.logjointlik(STATE)
}

compute.utilities.f.giv.y.f.t = function(STATE)
{
  compute.log.p.y.giv.f(STATE)
  compute.log.p.f.giv.theta(STATE)
}

compute.utilities.t.surrogate = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.update.S.theta.R.theta(STATE)
  compute.update.mtg(STATE)

  compute.update.f(STATE)
  compute.log.p.y.giv.f(STATE)
  compute.log.p.f.giv.theta(STATE)
  compute.log.p.g.giv.theta(STATE)
  compute.log.p.theta(STATE)
}

compute.logjointlik = function(STATE)
{
  STATE$logjointlik = STATE$log.p.y.giv.f + STATE$log.p.f.giv.theta + STATE$log.p.theta

  for(i in 1:SAMPLER$nblocks)
    STATE[[paste("LOGJOINTLIK.AR", i, sep="")]] = SAMPLER$blocks[[i]]$compute.logjointlik.AR(STATE)
}

compute.logjointlik.AR.t.surrogate = function(STATE)
{
  STATE$log.p.y.giv.f + STATE$log.p.g.giv.theta + STATE$log.p.theta
}

}


## **************************************************************************************************** PM
if(APPROACH == "PM") {

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

fun.log.p.y.giv.f = function(f)
  {
    res = sum(pnorm(DATA$y * f, log=T))
    
    if(is.nan(res)) res = -Inf
    res
  }

fun.log.p.f.giv.theta = function(f, STATE)
{
  nu = forwardsolve(STATE$L.chol.Q.mat, f) / sqrt(exp(STATE$psi.sigma))
  - sum(log(diag(STATE$L.chol.Q.mat))) - DATA$n / 2 * STATE$psi.sigma - 0.5 * crossprod(nu) - DATA$n/2 * log(2*pi)
}

compute.log.p.y.giv.theta.pseudo.marginal = function(STATE)
{
  if(APPROXIMATION == "EP") approx.model = EP.NOGRAD.MEAN.COV(STATE$psi.sigma, STATE$psi.tau, MAXIT=1)
  if(APPROXIMATION == "LA") approx.model = LA(STATE)

  f.pseudo = approx.model$L.Sigma %*% matrix(rnorm(DATA$n*NPSEUDO), DATA$n, NPSEUDO) + c(approx.model$mu)
  
  tmp = forwardsolve(approx.model$L.Sigma, f.pseudo - c(approx.model$mu))
  crossprodtmp = apply(tmp^2, 2, sum)
  log.pseudo.weights = apply(f.pseudo, 2, fun.log.p.y.giv.f) + apply(f.pseudo, 2, fun.log.p.f.giv.theta, STATE) + sum(log(diag(approx.model$L.Sigma))) + 0.5 * crossprodtmp + DATA$n/2 * log(2*pi) - log(NPSEUDO)
  
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

}

## ## **************************************************************************************************** 
## ## **************************************************************************************************** 



## **************************************************************************************************** PM
if(APPROACH == "PM.AIS") {

## compute.initial.utilities = function(STATE)
## {
##   compute.L.chol.Q.mat(STATE)

##   compute.log.p.y.giv.theta.pseudo.marginal(STATE)
##   compute.log.p.theta(STATE)
  
##   compute.logjointlik(STATE)
## }
 
compute.utilities.pseudo.marginal = function(STATE)
{
  compute.L.chol.Q.mat(STATE)
  compute.log.p.y.giv.theta.pseudo.marginal(STATE)
  compute.log.p.theta(STATE)
}

fun.log.p.y.giv.f = function(f)
  {
    res = sum(pnorm(DATA$y * f, log=T))
    
    if(is.nan(res)) res = -Inf
    res
  }

fun.log.p.f.giv.theta = function(f, STATE)
{
  nu = forwardsolve(STATE$L.chol.Q.mat, f) / sqrt(exp(STATE$psi.sigma))
  - sum(log(diag(STATE$L.chol.Q.mat))) - DATA$n / 2 * STATE$psi.sigma - 0.5 * crossprod(nu) - DATA$n/2 * log(2*pi)
}

compute.log.p.y.giv.theta.pseudo.marginal = function(STATE)
{
##   if(APPROXIMATION == "EP") approx.model = EP.NOGRAD.MEAN.COV(STATE$psi.sigma, STATE$psi.tau)
##   if(APPROXIMATION == "LA") approx.model = LA(STATE)

##   f.pseudo = approx.model$L.Sigma %*% matrix(rnorm(DATA$n*NPSEUDO), DATA$n, NPSEUDO) + c(approx.model$mu)

   ############ ANNEALING FROM THE PRIOR
  if(FROM.PRIOR == T) {
    log.pseudo.weights = rep(0, NPSEUDO)
    for(jjj in 1:NPSEUDO) {
      fpseudo = exp(STATE$psi.sigma/2) * STATE$L.chol.Q.mat %*% rnorm(DATA$n)
      log.pseudo.weights[jjj] = log.pseudo.weights[jjj] + (BETA.VECT[NBETA-1] - BETA.VECT[NBETA]) * fun.log.p.y.giv.f(fpseudo)

      if(NBETA > 2) {
        for(bbb in (NBETA-1):2) {
          for(rrr in 1:N.ELLSS.ANNEAL) fpseudo = ELLIPTICALSS.AIS(fpseudo, BETA.VECT[bbb], STATE)
          log.pseudo.weights[jjj] = log.pseudo.weights[jjj] + (BETA.VECT[bbb-1] - BETA.VECT[bbb]) * fun.log.p.y.giv.f(fpseudo)
        }
      }
    }

    res = LOGSUM.VECT(log.pseudo.weights) - log(NPSEUDO)
    if(is.nan(res)) res = -Inf
    
    STATE$log.p.y.giv.theta.pseudo = res
  }
  
  ################################## END

   ############ ANNEALING FROM THE LAPLACE APPROXIMATION
  if(FROM.PRIOR == F) {
  approx.model <<- LA(STATE)
  approx.model.normalization <<- - sum(log(diag(approx.model$L.Sigma))) - DATA$n/2*log(2*pi)
##   approx.model <<- EP.NOGRAD.MEAN.COV(STATE$psi.sigma, STATE$psi.tau)
  
  log.pseudo.weights = rep(0, NPSEUDO)
  for(jjj in 1:NPSEUDO) {
    fpseudo = approx.model$L.Sigma %*% rnorm(DATA$n) + approx.model$mu
    tmp = forwardsolve(approx.model$L.Sigma, fpseudo - c(approx.model$mu))
##     log.pseudo.weights[jjj] = log.pseudo.weights[jjj] + (BETA.VECT[NBETA-1] - BETA.VECT[NBETA]) * ( fun.log.p.y.giv.f(fpseudo) + fun.log.p.f.giv.theta(fpseudo, STATE) ) + (BETA.VECT[NBETA] - BETA.VECT[NBETA-1]) * (- sum(log(diag(approx.model$L.Sigma))) - 0.5 * crossprod(tmp) - DATA$n/2*log(2*pi) )
    log.pseudo.weights[jjj] = log.pseudo.weights[jjj] + (BETA.VECT[NBETA-1] - BETA.VECT[NBETA]) * ( fun.log.p.y.giv.f(fpseudo) + fun.log.p.f.giv.theta(fpseudo, STATE) ) + (BETA.VECT[NBETA] - BETA.VECT[NBETA-1]) * ( - 0.5 * crossprod(tmp) + approx.model.normalization )
    
    if(NBETA > 2) {
      for(bbb in (NBETA-1):2) {
        for(rrr in 1:N.ELLSS.ANNEAL) fpseudo = ELLIPTICALSS.AIS.LA(fpseudo, BETA.VECT[bbb], STATE)
        tmp = forwardsolve(approx.model$L.Sigma, fpseudo - c(approx.model$mu))
##         log.pseudo.weights[jjj] = log.pseudo.weights[jjj] + (BETA.VECT[bbb-1] - BETA.VECT[bbb]) * ( fun.log.p.y.giv.f(fpseudo) + fun.log.p.f.giv.theta(fpseudo, STATE) ) + (BETA.VECT[bbb] - BETA.VECT[bbb-1]) * ( - sum(log(diag(approx.model$L.Sigma))) - 0.5 * crossprod(tmp) - DATA$n/2*log(2*pi) )
        log.pseudo.weights[jjj] = log.pseudo.weights[jjj] + (BETA.VECT[bbb-1] - BETA.VECT[bbb]) * ( fun.log.p.y.giv.f(fpseudo) + fun.log.p.f.giv.theta(fpseudo, STATE) ) + (BETA.VECT[bbb] - BETA.VECT[bbb-1]) * ( - 0.5 * crossprod(tmp) + approx.model.normalization )
      }
    }
  }
  
  res = LOGSUM.VECT(log.pseudo.weights) - log(NPSEUDO)
  if(is.nan(res)) res = -Inf

  STATE$log.p.y.giv.theta.pseudo = res
}
  
##   tmp = forwardsolve(approx.model$L.Sigma, f.pseudo - c(approx.model$mu))
##   crossprodtmp = apply(tmp^2, 2, sum)
##   log.pseudo.weights = apply(f.pseudo, 2, fun.log.p.y.giv.f) + apply(f.pseudo, 2, fun.log.p.f.giv.theta, STATE) + sum(log(diag(approx.model$L.Sigma))) + 0.5 * crossprodtmp - log(NPSEUDO)
  
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




}
