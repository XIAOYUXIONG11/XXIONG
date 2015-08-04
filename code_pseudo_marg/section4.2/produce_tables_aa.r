## File to produce the rows of tables 1 and 3 (only for the AA and SURR approaches) - note that you first need to run the chains using run_experiment.r
## For the PM rows of tables 1 and 3 please run produce_tables.r


pdf.options(width=10,height=10,pointsize=32)
ps.options(width=10,height=10,paper="special",horizontal=F, pointsize=32)


library(coda)

NSEEDS = 10

## if APPROACH == "SURR" produce the rows of the table for the SURR method, otherwise produce the row of the table for the AA method
SCHEME.SAMPLER.HYPER = "MH"

for(NTRAINING in c(50, 200))
  {
    for(DTRAINING in c(2, 10))
      {
        for(APPROACH in c("AA", "SURR"))
          {
            for(COVARIANCE.MODE in c("H", "D"))
              {
                
                zz.acceptance = rep(0, NSEEDS)
                zz.ess = rep(0, NSEEDS)

                if(COVARIANCE.MODE == "H") tabellone = array(0, dim=c(15000, NSEEDS, 2))
                if(COVARIANCE.MODE == "D") tabellone = array(0, dim=c(15000, NSEEDS, DTRAINING + 1))

                for(SEED in 1:NSEEDS)
                  {
                    if(APPROACH == "AA")
                      {
                        ## Specify the directory where you have the results of your MCMC runs - AA method
                        FILESAVECHAINS.psi.tau = paste("aa/CHAINS_PSITAU_COV_", COVARIANCE.MODE, "_N_", NTRAINING, "_D_", DTRAINING, "_SEED_", SEED, ".txt", sep="")
                        FILESAVECHAINS.psi.sigma = paste("aa/CHAINS_PSISIGMA_COV_", COVARIANCE.MODE, "_N_", NTRAINING, "_D_", DTRAINING, "_SEED_", SEED, ".txt", sep="")
                      }
                    
                    if(APPROACH == "SURR")
                      {
                        ## Specify the directory where you have the results of your MCMC runs - SURR method
                        FILESAVECHAINS.psi.tau = paste("surr/CHAINS_PSITAU_COV_", COVARIANCE.MODE, "_N_", NTRAINING, "_D_", DTRAINING, "_SEED_", SEED, ".txt", sep="")
                        FILESAVECHAINS.psi.sigma = paste("surr/CHAINS_PSISIGMA_COV_", COVARIANCE.MODE, "_N_", NTRAINING, "_D_", DTRAINING, "_SEED_", SEED, ".txt", sep="")
                      }
                    
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
                
                if(APPROACH == "AA") cat("AA", " & $", pr1, "$ $(", pr2, ")$ & $", pr3, "$ & $", pr3.2, "$ & $", pr3.3, "$ & $", pr3.4, "$ ", sep="")
                if(APPROACH == "SURR") cat("SURR", " & $", pr1, "$ $(", pr2, ")$ & $", pr3, "$ & $", pr3.2, "$ & $", pr3.3, "$ & $", pr3.4, "$ ", sep="")
                cat(" & $", pr4, "$ $(", pr5, ")$", "  \\\\\n\n", sep="")
                
              }
            
  }
}
}
