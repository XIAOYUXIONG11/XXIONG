## Probit GP classification - Type II Maximum likelihood for optimizing the hyper-parameters.
## In the case of one length-scale parameters, the two hyper-parameters are optimized using grid search

APPROACH = "NONE"
source("../../functions/GP.FUNCTIONS.r")
source("../../functions/EP.r")

OMEGA = 1e-6

if(COVARIANCE.MODE == "H") NLS = 1
if(COVARIANCE.MODE == "D") NLS = D


PAR.VAL = rep(-Inf, NLS+1)

fun.to.optim = function(ppp)
  {
    if(prod(ppp == PAR.VAL))
      {
        return(MODEL.OPTIM$log.Z.ep)
      }
    
    assign("PAR.VAL", ppp, envir=.GlobalEnv)
    tmp = EP(psi.sigma=ppp[1], psi.tau=ppp[-1])
    assign("MODEL.OPTIM", tmp, envir=.GlobalEnv)
    
    tmp$log.Z.ep
  }

grad.fun.to.optim = function(ppp)
  {
    if(prod(ppp == PAR.VAL))
      {
        return(MODEL.OPTIM$grad.log.Z.ep)
      }
    
    assign("PAR.VAL", ppp, envir=.GlobalEnv)
    tmp = EP(psi.sigma=ppp[1], psi.tau=ppp[-1])
    assign("MODEL.OPTIM", tmp, envir=.GlobalEnv)
    
    tmp$grad.log.Z.ep
  }

if(COVARIANCE.MODE == "D")
  {
    best.model = list()
    best.model$value = Inf
    for(i in 1:10)
      {
        cat("Iteration", i, "   ")
        tmp = optim(runif(NLS+1, min=c(1,rep(-3, NLS)), max=c(3, rep(-1,NLS))), fun.to.optim, control=list(maxit=1000, fnscale=-1), method="L-BFGS-B", gr = grad.fun.to.optim)
        cat(tmp$value, "    par", tmp$par, "\n")
        ## tmp = optim(runif(NLS+1, min=c(1,rep(-3, NLS)), max=c(3, rep(-1,NLS))), fun.to.optim, control=list(maxit=1000, fnscale=-1), method="CG", gr = grad.fun.to.optim)
        if(tmp$value < best.model$value) best.model = tmp
      }
    
    cat("\n\n")
    
    ep.model = EP(psi.sigma=best.model$par[1], psi.tau=best.model$par[-1])
    predictions.test = predict.ep(best.model$par[1], best.model$par[-1], ep.model$nu.tilde.vect, ep.model$tau.tilde.vect)
  }

if(COVARIANCE.MODE == "H")
  {
    log.sigma.vect = seq(from=-1, to=12, length.out=15)
    log.tau.vect = seq(from=-4, to=6, length.out=15)
    
    best.value = -Inf
    best.par = NULL
    for(hhh in 1:length(log.sigma.vect))
      {
        for(jjj in 1:length(log.tau.vect))
          {
            tmp = fun.to.optim(c(log.sigma.vect[hhh], log.tau.vect[jjj]))
            cat(tmp, "\n")

            if(tmp > best.value)
              {
                best.value = tmp
                best.par = c(log.sigma.vect[hhh], log.tau.vect[jjj])
              }
          }
      }
    cat("\n")
    
    ep.model = EP(psi.sigma=best.par[1], psi.tau=best.par[-1])
    predictions.test = predict.ep(best.par[1], best.par[-1], ep.model$nu.tilde.vect, ep.model$tau.tilde.vect)

    print(best.par)
  }
