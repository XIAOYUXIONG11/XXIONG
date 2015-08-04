## **************************************************************************************************** 
## ************************************************** Laplace Approximation - Newton-Raphson
## **************************************************************************************************** 

## Putting the break before updating "converged" ensures that another iteration is performed before leaving the loop, so that W, B, L are updated and computed according to f.vect
LA = function(STATE)
  {
    f.vect = rep(0, DATA$n)

    converged = F
    for(iiii in 1:1000)
      {
        diag.W.mat = c(exp(2*dnorm(f.vect, log=T) - 2*pnorm(DATA$y * f.vect, log=T)) + DATA$y * f.vect * exp(dnorm(f.vect, log=T) - pnorm(DATA$y * f.vect, log=T)))

        B.mat = diag(DATA$n) + t(t(sqrt(diag.W.mat) * STATE$Q.mat * exp(STATE$psi.sigma)) * sqrt(diag.W.mat))
        L.chol = t(chol(B.mat))
        assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+1, envir=.GlobalEnv)
        if(converged == T) break

        grad.wrt.f = DATA$y * exp(dnorm(f.vect, log=T) - pnorm(DATA$y * f.vect, log=T))
        
        b.vect = diag.W.mat * f.vect + grad.wrt.f
        a.vect = b.vect - sqrt(diag.W.mat) * backsolve(t(L.chol), forwardsolve(L.chol, sqrt(diag.W.mat) * (exp(STATE$psi.sigma) * STATE$Q.mat %*% b.vect)))

        f.new = exp(STATE$psi.sigma) * STATE$Q.mat %*% a.vect
        
        if(mean((f.new - f.vect)^2) < 1e-4) converged = T 

        f.vect = f.new
      }

    res = list()
    
    res$mu = f.new
    
    tmp = forwardsolve(L.chol, (STATE$Q.mat) * (sqrt(diag.W.mat) * exp(STATE$psi.sigma)))
    SIGMA.F = exp(STATE$psi.sigma) * STATE$Q.mat - t(tmp) %*% tmp
    assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+2, envir=.GlobalEnv)
    res$L.Sigma = t(chol(SIGMA.F))
    assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+1, envir=.GlobalEnv)

    res
  }


LA.Z = function(STATE)
  {
    f.vect = rep(0, DATA$n)

    converged = F
    for(iiii in 1:1000)
      {
        diag.W.mat = c(exp(2*dnorm(f.vect, log=T) - 2*pnorm(DATA$y * f.vect, log=T)) + DATA$y * f.vect * exp(dnorm(f.vect, log=T) - pnorm(DATA$y * f.vect, log=T)))

        B.mat = diag(DATA$n) + t(t(sqrt(diag.W.mat) * STATE$Q.mat * exp(STATE$psi.sigma)) * sqrt(diag.W.mat))
        L.chol = t(chol(B.mat))
        assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+1, envir=.GlobalEnv)
        if(converged == T) break

        grad.wrt.f = DATA$y * exp(dnorm(f.vect, log=T) - pnorm(DATA$y * f.vect, log=T))
        
        b.vect = diag.W.mat * f.vect + grad.wrt.f
        a.vect = b.vect - sqrt(diag.W.mat) * backsolve(t(L.chol), forwardsolve(L.chol, sqrt(diag.W.mat) * (exp(STATE$psi.sigma) * STATE$Q.mat %*% b.vect)))

        f.new = exp(STATE$psi.sigma) * STATE$Q.mat %*% a.vect
        
        if(mean((f.new - f.vect)^2) < 1e-4) converged = T 

        f.vect = f.new
      }

    res = list()
    
    res$mu = f.new
    
    tmp = forwardsolve(L.chol, (STATE$Q.mat) * (sqrt(diag.W.mat) * exp(STATE$psi.sigma)))
    SIGMA.F = exp(STATE$psi.sigma) * STATE$Q.mat - t(tmp) %*% tmp
    assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+2, envir=.GlobalEnv)
    res$L.Sigma = t(chol(SIGMA.F))
    assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+1, envir=.GlobalEnv)

    tmp = forwardsolve(STATE$L.chol.Q.mat, f.new)
    res$log.Z.LA = -0.5 * crossprod(c(tmp)) + sum(pnorm(DATA$y * f.new, log=T)) - sum(log(diag(L.chol)))

    res
  }
