GENERAL.SAMPLER = function()
{
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
  
  
  ## Start the sampling procedure
  hhhh = 2
  for(ind.sampling in 2:SAMPLER$nsamples) {
    
    for(kkk in 1:SAMPLER$save.every) {
      
      ## Apply the transitions to the chain
      for(i in 1:SAMPLER$nblocks) {
        SAMPLES$BLOCK.INDEX = i
        if(APPROACH == "SURR") {
          if(i == 2) {
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

    hhhh = hhhh + 1
    
    if( hhhh > MAXTOSTORE ) {
      hhhh = 1
      
      if(CHECK.MODE == F)
        {
          if(APPROACH != "PM") write.table(formatC(SAMPLES$f, format="f", digits=3), file=FILESAVECHAINS.f, append=T, row.names=F, col.names=F, quote=F)
          write.table(formatC(SAMPLES$psi.sigma, format="f", digits=3), file=FILESAVECHAINS.psi.sigma, append=T, row.names=F, col.names=F, quote=F)
          write.table(formatC(SAMPLES$psi.tau, format="f", digits=3), file=FILESAVECHAINS.psi.tau, append=T, row.names=F, col.names=F, quote=F)
        }
    }
    
    ## Every nbatch steps print the acceptance rate
    if( (ind.sampling %% SAMPLER$nbatch ) == 0) {
      acceptance.rate = count.accepted / (SAMPLER$nbatch * SAMPLER$save.every )
      
      cat(ind.sampling/SAMPLER$nsamples*100, "%\t", sep="")
      cat(acceptance.rate, sep="\t")
      cat("\n")
      
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
          }
        }
      }
      count.accepted[] = 0
    }
    SAMPLES$IND.SAMPLING = ind.sampling
  }

  write.table(NUMBER.CALLS.N3, file=FILESAVE.NON3, append=T, row.names=F, col.names=F, quote=F)

  res = list()
  
  res
}
