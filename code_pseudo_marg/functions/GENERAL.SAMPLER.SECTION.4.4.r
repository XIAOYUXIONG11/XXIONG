GENERAL.SAMPLER = function()
{
  if(WITH.PREDICTIONS == T) assign("TEST.PREDICTIONS.MCMC", matrix(0, DATA$ntest, 0), envir=.GlobalEnv)
  
  SAMPLES <<- new.env()
  for(i in 1:length(MODEL$NAMES.GROUPS)) {
    SAMPLES[[MODEL$NAMES.GROUPS[i]]] = matrix(0, MAXTOSTORE, MODEL$SIZE.GROUPS[i])
    SAMPLES[[MODEL$NAMES.GROUPS[i]]][1,] = SAMPLER$init.par[[MODEL$NAMES.GROUPS[i]]]
  }

  cat("\n", "GENERAL SAMPLER\n", sep="")
  for(i in 1:SAMPLER$nblocks)
    cat("Sampling:", SAMPLER$blocks[[i]]$group, "using", SAMPLER$blocks[[i]]$method, "\n")

  reset.number.calls()

  ## Define chain as a structure containing the current state and variables related to the chain
  count.accepted = rep(0, SAMPLER$nblocks)
  plot.acceptance.rate = matrix(0, nrow=0, ncol=SAMPLER$nblocks)
  plot.logjointlikelihood = c()

  CURRENT.STATE <<- new.env()
  for(i in 1:length(MODEL$NAMES.GROUPS))
    CURRENT.STATE[[MODEL$NAMES.GROUPS[i]]] = SAMPLER$init.par[[MODEL$NAMES.GROUPS[i]]]
  PROPOSED.STATE <<- new.env()
  
  SAMPLER$compute.initial.utilities(CURRENT.STATE)

  if(APPROACH == "SURR") {
    SURROGATE.G <<- sqrt(CURRENT.STATE$Stheta) * rnorm(DATA$n) + CURRENT.STATE$f
    compute.update.mtg(CURRENT.STATE)
    compute.update.nu(CURRENT.STATE)
    
    compute.log.p.g.giv.theta(CURRENT.STATE)
    compute.logjointlik(CURRENT.STATE)
  }
  
  if(CHECK.MODE == T) est.mean.prior <<- exp(c(SAMPLES$psi.sigma[1], SAMPLES$psi.tau[1,]))
  if(CHECK.MODE == T) est.var.prior <<- rep(0, MODEL$NLS + 1)

  ## Start the sampling procedure
  hhhh = 2
  for(ind.sampling in 2:SAMPLER$nsamples) {
    plot.logjointlikelihood = c(plot.logjointlikelihood, CURRENT.STATE$logjointlik)
    
    for(kkk in 1:SAMPLER$save.every) {
      if(CHECK.MODE == T) {
        ## Generate the data, given the parameters
        if(LGM == "LOGREG") DATA$y = rbinom(DATA$n, 1, CURRENT.STATE$logistic.vect)
        if(LGM == "LOGCOX") DATA$y = rpois(DATA$n, exp(CURRENT.STATE$f))
        if(LGM == "COPULA") DATA$y = rnorm(DATA$n, 0, sd=exp(CURRENT.STATE$f))
        if(LGM == "ORDREG") {
          DATA$y = c()
          for(i in 1:DATA$n)
            {
              p.y = CURRENT.STATE$f
              TMP = diff(pnorm((B.INT - CURRENT.STATE$f[i]) / SIGMA.DELTA))
              DATA$y = c(DATA$y, sample(c(1:(length(B.INT)-1)), 1, replace=F, TMP))
            }
        }
        
        SAMPLER$compute.initial.utilities(CURRENT.STATE)
      }

      ## Apply the transitions to the chain
      for(i in 1:SAMPLER$nblocks) {
        SAMPLES$BLOCK.INDEX = i

        if(APPROACH == "SURR") {
          if("psi.tau" %in% BLOCKING[[i]]$group) {
            SURROGATE.G <<- sqrt(CURRENT.STATE$Stheta) * rnorm(DATA$n) + CURRENT.STATE$f
            compute.update.mtg(CURRENT.STATE)
            compute.update.nu(CURRENT.STATE)
            
            compute.log.p.g.giv.theta(CURRENT.STATE)
            compute.logjointlik(CURRENT.STATE)
          }
        }

        
        for(nested.repetitions in 1:SAMPLER$blocks[[i]]$repetitions)  {
          A = get(SAMPLER$blocks[[i]]$method)()
          u = -rexp(1, 1)

          if(u < A) {
            for(iii in ls(PROPOSED.STATE))
              CURRENT.STATE[[iii]] = PROPOSED.STATE[[iii]]
            
            count.accepted[i] = count.accepted[i]+1/SAMPLER$blocks[[i]]$repetitions
          }
        }
      }
    }

    for(i in 1:length(MODEL$NAMES.GROUPS))
      SAMPLES[[MODEL$NAMES.GROUPS[i]]][hhhh,] <<- CURRENT.STATE[[MODEL$NAMES.GROUPS[i]]]

    if(WITH.PREDICTIONS == T)
      {
        if(ind.sampling > SAMPLER$nburnin)
          {
            K.star.star = rep(0, DATA$ntest)
            K.star = exp(CURRENT.STATE$psi.sigma) * q.fun.xy(DATA$X, DATA$X.TEST, exp(CURRENT.STATE$psi.tau))
            for(i in 1:DATA$ntest) K.star.star[i] = exp(CURRENT.STATE$psi.sigma) * q.fun.xx(matrix(DATA$X.TEST[i,], nrow=1), exp(CURRENT.STATE$psi.sigma), exp(CURRENT.STATE$psi.tau))
            
            Lf = forwardsolve(exp(CURRENT.STATE$psi.sigma/2) * CURRENT.STATE$L.chol.Q.mat, CURRENT.STATE$f)
            LKstar = forwardsolve(CURRENT.STATE$L.chol.Q.mat, K.star)

            mstar = t(LKstar) %*% Lf ## t(K.star) %*% solve(K) %*% f
            vstar = K.star.star - apply(LKstar^2, 2, sum) ## K.star.star - apply((solve(K) %*% K.star) * K.star, 2, sum)
            vstar[vstar < 0] = 1e-9
    
            fstar = matrix(rnorm(DATA$ntest*100, mean=mstar, sd=sqrt(vstar)), nrow=DATA$ntest)
            
            p.y.giv.f.star = apply(exp(pnorm(fstar, log=T)), 1, mean)

            TEST.PREDICTIONS.MCMC <<- cbind(TEST.PREDICTIONS.MCMC, p.y.giv.f.star)
          }
      }

    if(CHECK.MODE == T) {
      est.var.prior <<- (est.var.prior * (ind.sampling-1) / (ind.sampling) + (exp(c(SAMPLES$psi.sigma[hhhh], SAMPLES$psi.tau[hhhh,])) - est.mean.prior)^2 * (ind.sampling-1) / (ind.sampling)^2)
      est.mean.prior <<- est.mean.prior * (ind.sampling-1) / (ind.sampling) + exp(c(SAMPLES$psi.sigma[hhhh], SAMPLES$psi.tau[hhhh,])) / (ind.sampling)
    }

    hhhh = hhhh + 1
    
    if( hhhh > MAXTOSTORE ) {
      hhhh = 1
      
##       if(CHECK.MODE == F)
##         {
##       write.table(formatC(SAMPLES$f, format="f", digits=3), file=FILESAVECHAINS.f, append=T, row.names=F, col.names=F, quote=F)
          write.table(formatC(SAMPLES$psi.sigma, format="f", digits=3), file=FILESAVECHAINS.psi.sigma, append=T, row.names=F, col.names=F, quote=F)
          write.table(formatC(SAMPLES$psi.tau, format="f", digits=3), file=FILESAVECHAINS.psi.tau, append=T, row.names=F, col.names=F, quote=F)
##         }
    }
    
    ## Every nbatch steps print the acceptance rate
    if( (ind.sampling %% SAMPLER$nbatch ) == 0) {
##         if(ind.sampling > (SAMPLER$nburnin+1000)) {
##           tmp = bmmat(t(TEST.PREDICTIONS.MCMC))
##           row.names(tmp) = NULL
          
##           predictions.test = tmp[,1]
##           predictions.test.se = tmp[,2]
          
##           if(max(predictions.test.se) < 0.05) {
##             write.table(formatC(SAMPLES$f[1:ind.sampling,], format="f", digits=3), file=FILESAVECHAINS.f, append=T, row.names=F, col.names=F, quote=F)
##             write.table(formatC(SAMPLES$psi.sigma[1:ind.sampling,], format="f", digits=3), file=FILESAVECHAINS.psi.sigma, append=T, row.names=F, col.names=F, quote=F)
##             write.table(formatC(SAMPLES$psi.tau[1:ind.sampling,], format="f", digits=3), file=FILESAVECHAINS.psi.tau, append=T, row.names=F, col.names=F, quote=F)
            
##             assign("predictions.test", predictions.test, envir=.GlobalEnv)
##             assign("predictions.test.se", predictions.test.se, envir=.GlobalEnv)
##             assign("index.mcmc", ind.sampling, envir=.GlobalEnv)
##             return(0)
##           }
##         }
      
      acceptance.rate = count.accepted / (SAMPLER$nbatch * SAMPLER$save.every )
      plot.acceptance.rate = rbind(plot.acceptance.rate, acceptance.rate)
      
      cat(ind.sampling/SAMPLER$nsamples*100, "%\t", sep="")
      cat(acceptance.rate, sep="\t")
      cat("\n")
      
      ## Interrupt if the check is passed
      if(CHECK.MODE == T) {
          if( ( max(abs(est.mean.prior - true.mean.prior)) < 1e-1 ) & ( max(abs(sqrt(est.var.prior) - sqrt(true.var.prior))) < 1e-1) ) return(0)
        }

      if(ADAPTIVE == T) {
        if(ind.sampling < SAMPLER$nburnin) {
          ADAPTIVE.FACTOR = 1.1
          ADAPTIVE.FACTOR = 5 - min(ind.sampling, SAMPLER$nburnin/5) / (SAMPLER$nburnin/5) * 3.9
          print(ADAPTIVE.FACTOR)
          for(i in 1:SAMPLER$nblocks) {
            if(SAMPLER$blocks[[i]]$method == "MH") {
              if (acceptance.rate[i] < 0.15) SAMPLER$blocks[[i]]$cov.proposal = SAMPLER$blocks[[i]]$cov.proposal / ADAPTIVE.FACTOR
              if (acceptance.rate[i] > 0.30) SAMPLER$blocks[[i]]$cov.proposal = SAMPLER$blocks[[i]]$cov.proposal * ADAPTIVE.FACTOR
              SAMPLER$blocks[[i]]$L.chol.cov.proposal = t(chol(SAMPLER$blocks[[i]]$cov.proposal))
            }
            if(SAMPLER$blocks[[i]]$method == "HMC" || SAMPLER$blocks[[i]]$method == "HMCinv.f") {
              if (acceptance.rate[i] < 0.7) SAMPLER$blocks[[i]]$max.distance.travel = SAMPLER$blocks[[i]]$max.distance.travel / ADAPTIVE.FACTOR
              if (acceptance.rate[i] > 0.9) SAMPLER$blocks[[i]]$max.distance.travel = SAMPLER$blocks[[i]]$max.distance.travel * ADAPTIVE.FACTOR
              SAMPLER$blocks[[i]]$epsilon = SAMPLER$blocks[[i]]$max.distance.travel / SAMPLER$blocks[[i]]$L.max
            }
            if(SAMPLER$blocks[[i]]$method == "SMMALA") {
              if (acceptance.rate[i] < 0.5) SAMPLER$blocks[[i]]$epsilon = SAMPLER$blocks[[i]]$epsilon / ADAPTIVE.FACTOR
              if (acceptance.rate[i] > 0.7) SAMPLER$blocks[[i]]$epsilon = SAMPLER$blocks[[i]]$epsilon * ADAPTIVE.FACTOR
            }
          }
        }
      }
      count.accepted[] = 0
    }
    SAMPLES$IND.SAMPLING = ind.sampling
  }
  
  res = list()
  res$ess = list()
##   for(i in 1:length(MODEL$NAMES.GROUPS))
##     res$ess[[MODEL$NAMES.GROUPS[i]]] = ess(matrix(SAMPLES[[MODEL$NAMES.GROUPS[i]]][(SAMPLER$nsamples/2):SAMPLER$nsamples,], ncol=dim(SAMPLES[[MODEL$NAMES.GROUPS[i]]])[2]))
  res$plot.acceptance.rate = plot.acceptance.rate
  res$plot.logjointlikelihood = plot.logjointlikelihood
  
##   res$number.calls = list()
##   res$number.calls$logjointlikelihood = NUMBER.CALLS.LOGJOINTLIKELIHOOD
##   res$number.calls$grad.logjointlikelihood = NUMBER.CALLS.GRAD.LOGJOINTLIKELIHOOD
##   res$number.calls$metric.tensor = NUMBER.CALLS.METRIC.TENSOR
##   res$number.calls$covariance = NUMBER.CALLS.COVARIANCE
##   res$number.calls$covariance.cholesky = NUMBER.CALLS.COVARIANCE.CHOLESKY

  res
}
