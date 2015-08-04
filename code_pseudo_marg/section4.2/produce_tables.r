## File to produce the rows of tables 1 and 3 (only for the PM approach) - note that you first need to run the chains using run_experiment.r
## For the AA and SURR rows of tables 1 and 3 please run produce_tables_aa.r

pdf.options(width=10,height=10,pointsize=32)
ps.options(width=10,height=10,paper="special",horizontal=F, pointsize=32)


## Specify here where the chains for the PM approach are
directory.results = "mcmc_pseudo/"

library(coda)

NSEEDS = 10

SCHEME.SAMPLER.HYPER = "MH"

for(NTRAINING in c(50, 200))
  {
    for(DTRAINING in c(2, 10))
      {
        for(APPROXIMATION in c("LA", "EP"))
          {
            for(NPSEUDO in c(1, 16, 64))
              {
                for(COVARIANCE.MODE in c("H", "D"))
                  {
                zz.acceptance = rep(0, NSEEDS)
                zz.ess = rep(0, NSEEDS)
                ##                 zz.non3 = rep(0, NSEEDS)
                if(COVARIANCE.MODE == "H") tabellone = array(0, dim=c(15000, NSEEDS, 2))
                if(COVARIANCE.MODE == "D") tabellone = array(0, dim=c(15000, NSEEDS, DTRAINING + 1))
                
                for(SEED in 1:NSEEDS)
                  {
                    FILESAVECHAINS.psi.tau = paste(directory.results, "CHAINS_PSITAU_", NTRAINING, "_", DTRAINING, "_", APPROXIMATION, "_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
                    FILESAVECHAINS.psi.sigma = paste(directory.results, "CHAINS_PSISIGMA_", NTRAINING, "_", DTRAINING, "_", APPROXIMATION, "_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
                    
                    tabellone[,SEED,1] = as.matrix(read.table(FILESAVECHAINS.psi.sigma))
                    tabellone[,SEED,-1] = as.matrix(read.table(FILESAVECHAINS.psi.tau))

                    aa = min(effectiveSize(tabellone[5001:15000,SEED,]))
                    zz.ess[SEED] = aa

                    aa = sum(diff(tabellone[5001:15000,SEED,][,1]) != 0)
                    zz.acceptance[SEED] = aa
                  }
                
                tmp = mcmc.list()
                for(i in 1:NSEEDS) tmp[[i]] = mcmc(tabellone[1:15000,i,])

                zz = gelman.plot(tmp, autoburnin=F, max.bin=299)
                ind.rhat = which(zz$last.iter %in% c(6001, 7001, 10001, 15000))
                RHATtmp = apply(zz$shrink[,,1], 1, max)[ind.rhat]
                                
                pr1 = formatC(mean(zz.ess), format="f", digits=0)
                pr2 = formatC(sd(zz.ess), format="f", digits=0)
                pr3 = formatC(RHATtmp[1], format="f", digits=2)
                pr3.2 = formatC(RHATtmp[2], format="f", digits=2)
                pr3.3 = formatC(RHATtmp[3], format="f", digits=2)
                pr3.4 = formatC(RHATtmp[4], format="f", digits=2)
                pr4 = formatC(mean(zz.acceptance/10000*100), format="f", digits=1)
                pr5 = formatC(sd(zz.acceptance/10000*100), format="f", digits=1)

                print(FILESAVECHAINS.psi.sigma)
                
                cat("PM ", APPROXIMATION, " $(", NPSEUDO, ")$", " & $", pr1, "$ $(", pr2, ")$ & $", pr3, "$ & $", pr3.2, "$ & $", pr3.3, "$ & $", pr3.4, "$ ", sep="")
                cat(" & $", pr4, "$ $(", pr5, ")$", "  \\\\\n\n", sep="")

                
              }
          }
            cat("\n\n", sep="")
      }
  }
}
