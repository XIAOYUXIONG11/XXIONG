## This file contains the EP algorithm for probit GP classification - the implementation follows the book Rasmussen and Williams (2006)
## The file contains a few different versions depending on what they return. All versions return nu.tilde.vect and tau.tilde.vect as well as:
## EP - returns log evidence and gradients
## EP.NOGRAD - returns the log evidence
## EP.NOGRAD.MEAN.COV - returns the log evidence and the mean and covariance of the Gaussian approximation to the posterior of the latent variables

## The function predict.ep returns the prediction on test data according to the approximation carried out by one of the EP functions

EP = function(psi.sigma, psi.tau, MAXIT=1000)
  {
    K = Sigma = exp(psi.sigma) * q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))

    nu.tilde.vect = rep(0, DATA$n)
    tau.tilde.vect = rep(0, DATA$n)

    mu.cavity.vect = rep(0, DATA$n)
    sigma2.cavity.vect = rep(0, DATA$n)
    
    mu.vect = rep(0, DATA$n)
    log.Z.ep = -Inf
    
    for(iteration in 1:MAXIT)
      {
        for(i in 1:DATA$n)
          {
            tau.cavity = 1/Sigma[i,i] - tau.tilde.vect[i]
            nu.cavity = mu.vect[i] / Sigma[i,i] - nu.tilde.vect[i]
            mu.cavity.vect[i] = nu.cavity / tau.cavity
            sigma2.cavity.vect[i] = 1/tau.cavity

            zi = DATA$y[i] * nu.cavity / tau.cavity / sqrt(1 + 1/tau.cavity)
            Z.hat.vect.i = pnorm(zi)
            mu.hat.vect.i = nu.cavity / tau.cavity + DATA$y[i] * dnorm(zi) / tau.cavity / pnorm(zi) / sqrt(1 + 1/tau.cavity)
            sigma2.hat.vect.i = 1 / tau.cavity - 1 / (tau.cavity)^2 * dnorm(zi) / (1 + 1/tau.cavity) / pnorm(zi) * (zi + dnorm(zi) / pnorm(zi))
            
            delta.tau.tilde = 1/sigma2.hat.vect.i - tau.cavity - tau.tilde.vect[i]
            tau.tilde.vect[i] = max(0, tau.tilde.vect[i] + delta.tau.tilde)
            nu.tilde.vect[i] = mu.hat.vect.i / sigma2.hat.vect.i - nu.cavity

            s.vect = Sigma[,i]
            Sigma = Sigma - 1/(1/delta.tau.tilde + Sigma[i,i]) * s.vect %*% t(s.vect)
            mu.vect = Sigma %*% nu.tilde.vect
          }
        L = t(chol(diag(DATA$n) + t(sqrt(tau.tilde.vect) * K) * sqrt(tau.tilde.vect)))   
        
        V = forwardsolve(L, sqrt(tau.tilde.vect) * K)
        Sigma = K - t(V) %*% V
        mu.vect = Sigma %*% nu.tilde.vect

        p1p4 = 0.5 * sum(log(1 + tau.tilde.vect*sigma2.cavity.vect)) - sum(log(diag(L)))
        
        Knu = K %*% nu.tilde.vect
        s2t = sigma2.cavity.vect * tau.tilde.vect
        p2p5 = - 0.5 * crossprod(forwardsolve(L, Knu * sqrt(tau.tilde.vect))) + 0.5 * crossprod(mu.cavity.vect, (tau.tilde.vect * mu.cavity.vect - 2 * nu.tilde.vect)/(1 + s2t)) - 0.5 * crossprod(nu.tilde.vect * sigma2.cavity.vect / (1 + s2t), nu.tilde.vect) + 0.5 * crossprod(Knu, nu.tilde.vect)

        p3 = sum(pnorm(DATA$y * mu.cavity.vect / sqrt(1 + sigma2.cavity.vect), log=T))
        
        new.log.Z.ep = p1p4 + p3 + p2p5

        if(abs(new.log.Z.ep - log.Z.ep) < 1e-3) break
        log.Z.ep = new.log.Z.ep
      }

    b.vect = nu.tilde.vect - backsolve(t(L), forwardsolve(L, Knu * sqrt(tau.tilde.vect))) * sqrt(tau.tilde.vect)
    R = b.vect %*% t(b.vect) - crossprod(forwardsolve(L, diag(sqrt(tau.tilde.vect))))

    grad.log.Z.ep = rep(0, NLS + 1)
    grad.log.Z.ep[1] = 0.5 * sum(diag(R %*% (K - diag(OMEGA, DATA$n))))
    for(i in 1:NLS)
      {
        C = exp(- 2 * psi.tau[i]) * (K - diag(rep(OMEGA, DATA$n))) * abs(outer(DATA$X[,i], DATA$X[,i], "-"))^2
        grad.log.Z.ep[i+1] = 0.5 * sum(diag(R %*% C)) ## Can optimize this!
      }

    list(nu.tilde.vect = nu.tilde.vect, tau.tilde.vect = tau.tilde.vect, log.Z.ep = new.log.Z.ep, grad.log.Z.ep=grad.log.Z.ep)
  }

EP.NOGRAD = function(psi.sigma, psi.tau, MAXIT=1000)
  {
    K = Sigma = exp(psi.sigma) * q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))

    nu.tilde.vect = rep(0, DATA$n)
    tau.tilde.vect = rep(0, DATA$n)

    mu.cavity.vect = rep(0, DATA$n)
    sigma2.cavity.vect = rep(0, DATA$n)
    
    mu.vect = rep(0, DATA$n)
    log.Z.ep = -Inf
    mu.old = mu.vect
    
    for(iteration in 1:MAXIT)
      {
        for(i in 1:DATA$n)
          {
            tau.cavity = 1/Sigma[i,i] - tau.tilde.vect[i]
            nu.cavity = mu.vect[i] / Sigma[i,i] - nu.tilde.vect[i]
            mu.cavity.vect[i] = nu.cavity / tau.cavity
            sigma2.cavity.vect[i] = 1/tau.cavity

            zi = DATA$y[i] * nu.cavity / tau.cavity / sqrt(1 + 1/tau.cavity)
            Z.hat.vect.i = pnorm(zi)
            mu.hat.vect.i = nu.cavity / tau.cavity + DATA$y[i] * dnorm(zi) / tau.cavity / pnorm(zi) / sqrt(1 + 1/tau.cavity)
            sigma2.hat.vect.i = 1 / tau.cavity - 1 / (tau.cavity)^2 * dnorm(zi) / (1 + 1/tau.cavity) / pnorm(zi) * (zi + dnorm(zi) / pnorm(zi))
            
            delta.tau.tilde = 1/sigma2.hat.vect.i - tau.cavity - tau.tilde.vect[i]
            tau.tilde.vect[i] = max(0, tau.tilde.vect[i] + delta.tau.tilde)
            nu.tilde.vect[i] = mu.hat.vect.i / sigma2.hat.vect.i - nu.cavity

            s.vect = Sigma[,i]
            Sigma = Sigma - 1/(1/delta.tau.tilde + Sigma[i,i]) * s.vect %*% t(s.vect)
            mu.vect = Sigma %*% nu.tilde.vect
          }
        L = t(chol(diag(DATA$n) + t(sqrt(tau.tilde.vect) * K) * sqrt(tau.tilde.vect)))   
        
        V = forwardsolve(L, sqrt(tau.tilde.vect) * K)
        Sigma = K - t(V) %*% V
        mu.vect = Sigma %*% nu.tilde.vect

        p1p4 = 0.5 * sum(log(1 + tau.tilde.vect*sigma2.cavity.vect)) - sum(log(diag(L)))
        
        Knu = K %*% nu.tilde.vect
        s2t = sigma2.cavity.vect * tau.tilde.vect
        p2p5 = - 0.5 * crossprod(forwardsolve(L, Knu * sqrt(tau.tilde.vect))) + 0.5 * crossprod(mu.cavity.vect, (tau.tilde.vect * mu.cavity.vect - 2 * nu.tilde.vect)/(1 + s2t)) - 0.5 * crossprod(nu.tilde.vect * sigma2.cavity.vect / (1 + s2t), nu.tilde.vect) + 0.5 * crossprod(Knu, nu.tilde.vect)

        p3 = sum(pnorm(DATA$y * mu.cavity.vect / sqrt(1 + sigma2.cavity.vect), log=T))
        
        new.log.Z.ep = p1p4 + p3 + p2p5
        
        if(mean((mu.vect - mu.old)^2) < 1e-4) break
        log.Z.ep = new.log.Z.ep
        mu.old = mu.vect
      }

    list(nu.tilde.vect = nu.tilde.vect, tau.tilde.vect = tau.tilde.vect, log.Z.ep = new.log.Z.ep)
  }


EP.NOGRAD.CONVERGENCE.ON.Z = function(psi.sigma, psi.tau, MAXIT=1000)
  {
    K = Sigma = exp(psi.sigma) * q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))

    nu.tilde.vect = rep(0, DATA$n)
    tau.tilde.vect = rep(0, DATA$n)

    mu.cavity.vect = rep(0, DATA$n)
    sigma2.cavity.vect = rep(0, DATA$n)
    
    mu.vect = rep(0, DATA$n)
    log.Z.ep = -Inf
    
    for(iteration in 1:MAXIT)
      {
        for(i in 1:DATA$n)
          {
            tau.cavity = 1/Sigma[i,i] - tau.tilde.vect[i]
            nu.cavity = mu.vect[i] / Sigma[i,i] - nu.tilde.vect[i]
            mu.cavity.vect[i] = nu.cavity / tau.cavity
            sigma2.cavity.vect[i] = 1/tau.cavity

            zi = DATA$y[i] * nu.cavity / tau.cavity / sqrt(1 + 1/tau.cavity)
            Z.hat.vect.i = pnorm(zi)
            mu.hat.vect.i = nu.cavity / tau.cavity + DATA$y[i] * dnorm(zi) / tau.cavity / pnorm(zi) / sqrt(1 + 1/tau.cavity)
            sigma2.hat.vect.i = 1 / tau.cavity - 1 / (tau.cavity)^2 * dnorm(zi) / (1 + 1/tau.cavity) / pnorm(zi) * (zi + dnorm(zi) / pnorm(zi))
            
            delta.tau.tilde = 1/sigma2.hat.vect.i - tau.cavity - tau.tilde.vect[i]
            tau.tilde.vect[i] = max(0, tau.tilde.vect[i] + delta.tau.tilde)
            nu.tilde.vect[i] = mu.hat.vect.i / sigma2.hat.vect.i - nu.cavity

            s.vect = Sigma[,i]
            Sigma = Sigma - 1/(1/delta.tau.tilde + Sigma[i,i]) * s.vect %*% t(s.vect)
            mu.vect = Sigma %*% nu.tilde.vect
          }
        L = t(chol(diag(DATA$n) + t(sqrt(tau.tilde.vect) * K) * sqrt(tau.tilde.vect)))   
        
        V = forwardsolve(L, sqrt(tau.tilde.vect) * K)
        Sigma = K - t(V) %*% V
        mu.vect = Sigma %*% nu.tilde.vect

        p1p4 = 0.5 * sum(log(1 + tau.tilde.vect*sigma2.cavity.vect)) - sum(log(diag(L)))
        
        Knu = K %*% nu.tilde.vect
        s2t = sigma2.cavity.vect * tau.tilde.vect
        p2p5 = - 0.5 * crossprod(forwardsolve(L, Knu * sqrt(tau.tilde.vect))) + 0.5 * crossprod(mu.cavity.vect, (tau.tilde.vect * mu.cavity.vect - 2 * nu.tilde.vect)/(1 + s2t)) - 0.5 * crossprod(nu.tilde.vect * sigma2.cavity.vect / (1 + s2t), nu.tilde.vect) + 0.5 * crossprod(Knu, nu.tilde.vect)

        p3 = sum(pnorm(DATA$y * mu.cavity.vect / sqrt(1 + sigma2.cavity.vect), log=T))
        
        new.log.Z.ep = p1p4 + p3 + p2p5

        if(abs(new.log.Z.ep - log.Z.ep) < 1e-3) break
        log.Z.ep = new.log.Z.ep
      }

    list(nu.tilde.vect = nu.tilde.vect, tau.tilde.vect = tau.tilde.vect, log.Z.ep = new.log.Z.ep)
  }


EP.NOGRAD.MEAN.COV = function(psi.sigma, psi.tau, MAXIT=1000)
  {
    K = Sigma = exp(psi.sigma) * q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))
    
    nu.tilde.vect = rep(0, DATA$n)
    tau.tilde.vect = rep(0, DATA$n)

    mu.cavity.vect = rep(0, DATA$n)
    sigma2.cavity.vect = rep(0, DATA$n)
    
    mu.vect = rep(0, DATA$n)
    log.Z.ep = -Inf
    mu.old = mu.vect
    
    for(iteration in 1:MAXIT)
      {
        for(i in 1:DATA$n)
          {
            tau.cavity = 1/Sigma[i,i] - tau.tilde.vect[i]
            nu.cavity = mu.vect[i] / Sigma[i,i] - nu.tilde.vect[i]
            mu.cavity.vect[i] = nu.cavity / tau.cavity
            sigma2.cavity.vect[i] = 1/tau.cavity

            zi = DATA$y[i] * nu.cavity / tau.cavity / sqrt(1 + 1/tau.cavity)
            Z.hat.vect.i = pnorm(zi)
            mu.hat.vect.i = nu.cavity / tau.cavity + DATA$y[i] * dnorm(zi) / tau.cavity / pnorm(zi) / sqrt(1 + 1/tau.cavity)
            sigma2.hat.vect.i = 1 / tau.cavity - 1 / (tau.cavity)^2 * dnorm(zi) / (1 + 1/tau.cavity) / pnorm(zi) * (zi + dnorm(zi) / pnorm(zi))
            
            delta.tau.tilde = 1/sigma2.hat.vect.i - tau.cavity - tau.tilde.vect[i]
            tau.tilde.vect[i] = max(0, tau.tilde.vect[i] + delta.tau.tilde)
            nu.tilde.vect[i] = mu.hat.vect.i / sigma2.hat.vect.i - nu.cavity

            s.vect = Sigma[,i]
            Sigma = Sigma - 1/(1/delta.tau.tilde + Sigma[i,i]) * s.vect %*% t(s.vect)
            mu.vect = Sigma %*% nu.tilde.vect
          }
        L = t(chol(diag(DATA$n) + t(sqrt(tau.tilde.vect) * K) * sqrt(tau.tilde.vect)))   
        V = forwardsolve(L, sqrt(tau.tilde.vect) * K)
        Sigma = K - t(V) %*% V

        assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+3, envir=.GlobalEnv) 
        
        mu.vect = Sigma %*% nu.tilde.vect

        p1p4 = 0.5 * sum(log(1 + tau.tilde.vect*sigma2.cavity.vect)) - sum(log(diag(L)))
        
        Knu = K %*% nu.tilde.vect
        s2t = sigma2.cavity.vect * tau.tilde.vect
        p2p5 = - 0.5 * crossprod(forwardsolve(L, Knu * sqrt(tau.tilde.vect))) + 0.5 * crossprod(mu.cavity.vect, (tau.tilde.vect * mu.cavity.vect - 2 * nu.tilde.vect)/(1 + s2t)) - 0.5 * crossprod(nu.tilde.vect * sigma2.cavity.vect / (1 + s2t), nu.tilde.vect) + 0.5 * crossprod(Knu, nu.tilde.vect)

        p3 = sum(pnorm(DATA$y * mu.cavity.vect / sqrt(1 + sigma2.cavity.vect), log=T))
        
        new.log.Z.ep = p1p4 + p3 + p2p5

        if(mean((mu.vect - mu.old)^2) < 1e-4) break
        log.Z.ep = new.log.Z.ep
        mu.old = mu.vect
      }

    assign("NUMBER.CALLS.N3", NUMBER.CALLS.N3+1, envir=.GlobalEnv) 
    
    list(nu.tilde.vect = nu.tilde.vect, tau.tilde.vect = tau.tilde.vect, log.Z.ep = new.log.Z.ep, mu = Sigma %*% nu.tilde.vect, L.Sigma = t(chol(Sigma)))
  }


EP.NOGRAD.MEAN.COV.CONVERGENCE.ON.Z = function(psi.sigma, psi.tau, MAXIT=1000)
  {
    K = Sigma = exp(psi.sigma) * q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))

    nu.tilde.vect = rep(0, DATA$n)
    tau.tilde.vect = rep(0, DATA$n)

    mu.cavity.vect = rep(0, DATA$n)
    sigma2.cavity.vect = rep(0, DATA$n)
    
    mu.vect = rep(0, DATA$n)
    log.Z.ep = -Inf
    
    for(iteration in 1:MAXIT)
      {
        for(i in 1:DATA$n)
          {
            tau.cavity = 1/Sigma[i,i] - tau.tilde.vect[i]
            nu.cavity = mu.vect[i] / Sigma[i,i] - nu.tilde.vect[i]
            mu.cavity.vect[i] = nu.cavity / tau.cavity
            sigma2.cavity.vect[i] = 1/tau.cavity

            zi = DATA$y[i] * nu.cavity / tau.cavity / sqrt(1 + 1/tau.cavity)
            Z.hat.vect.i = pnorm(zi)
            mu.hat.vect.i = nu.cavity / tau.cavity + DATA$y[i] * dnorm(zi) / tau.cavity / pnorm(zi) / sqrt(1 + 1/tau.cavity)
            sigma2.hat.vect.i = 1 / tau.cavity - 1 / (tau.cavity)^2 * dnorm(zi) / (1 + 1/tau.cavity) / pnorm(zi) * (zi + dnorm(zi) / pnorm(zi))
            
            delta.tau.tilde = 1/sigma2.hat.vect.i - tau.cavity - tau.tilde.vect[i]
            tau.tilde.vect[i] = max(0, tau.tilde.vect[i] + delta.tau.tilde)
            nu.tilde.vect[i] = mu.hat.vect.i / sigma2.hat.vect.i - nu.cavity

            s.vect = Sigma[,i]
            Sigma = Sigma - 1/(1/delta.tau.tilde + Sigma[i,i]) * s.vect %*% t(s.vect)
            mu.vect = Sigma %*% nu.tilde.vect
          }
        L = t(chol(diag(DATA$n) + t(sqrt(tau.tilde.vect) * K) * sqrt(tau.tilde.vect)))   
        
        V = forwardsolve(L, sqrt(tau.tilde.vect) * K)
        Sigma = K - t(V) %*% V
        mu.vect = Sigma %*% nu.tilde.vect

        p1p4 = 0.5 * sum(log(1 + tau.tilde.vect*sigma2.cavity.vect)) - sum(log(diag(L)))
        
        Knu = K %*% nu.tilde.vect
        s2t = sigma2.cavity.vect * tau.tilde.vect
        p2p5 = - 0.5 * crossprod(forwardsolve(L, Knu * sqrt(tau.tilde.vect))) + 0.5 * crossprod(mu.cavity.vect, (tau.tilde.vect * mu.cavity.vect - 2 * nu.tilde.vect)/(1 + s2t)) - 0.5 * crossprod(nu.tilde.vect * sigma2.cavity.vect / (1 + s2t), nu.tilde.vect) + 0.5 * crossprod(Knu, nu.tilde.vect)

        p3 = sum(pnorm(DATA$y * mu.cavity.vect / sqrt(1 + sigma2.cavity.vect), log=T))
        
        new.log.Z.ep = p1p4 + p3 + p2p5

        if(abs(new.log.Z.ep - log.Z.ep) < 1e-3) break
        log.Z.ep = new.log.Z.ep
      }

    list(nu.tilde.vect = nu.tilde.vect, tau.tilde.vect = tau.tilde.vect, log.Z.ep = new.log.Z.ep, mu = Sigma %*% nu.tilde.vect, L.Sigma = t(chol(Sigma)))
  }


predict.ep = function(psi.sigma, psi.tau, nu.tilde.vect, tau.tilde.vect)
  {
    K = exp(psi.sigma) * q.fun.xx(DATA$X, exp(psi.sigma), exp(psi.tau))
    K.star = exp(psi.sigma) * q.fun.xy(DATA$X, DATA$X.TEST, exp(psi.tau))

    L = t(chol(diag(DATA$n) + t(sqrt(tau.tilde.vect) * K) * sqrt(tau.tilde.vect)))

    E.fstar = t(K.star) %*% (nu.tilde.vect - backsolve(t(L), forwardsolve(L, ((K %*% nu.tilde.vect) * sqrt(tau.tilde.vect)))) * sqrt(tau.tilde.vect))

    V.fstar = rep(0, DATA$ntest)
    for(i in 1:DATA$ntest)
      V.fstar[i] = exp(psi.sigma) * q.fun.xx(matrix(DATA$X.TEST[i,], nrow=1), exp(psi.sigma), exp(psi.tau)) - crossprod(forwardsolve(L, K.star[,i] * sqrt(tau.tilde.vect)))

    pnorm(E.fstar / sqrt(1 + V.fstar))
  }
