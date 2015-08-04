## File to produce the rows of tables 2 (only for the PM approach) - note that you first need to run the chains using run_experiment.r

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
                    zz.non3 = rep(0, NSEEDS)
                    
                    for(SEED in 1:NSEEDS)
                      {
                        FILESAVE.NON3 = paste(directory.results, "NON3_", NTRAINING, "_", DTRAINING, "_", APPROXIMATION, "_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_HYP_", SCHEME.SAMPLER.HYPER, "_SEED_", SEED, ".txt", sep="")
                    
                        aa = read.table(FILESAVE.NON3)[[1]]
                        zz.non3[SEED] = aa
                      }
                
                    pr1 = formatC(mean(zz.non3) / 15000, format="f", digits=1)
                    pr2 = formatC(sd(zz.non3) / 15000, format="f", digits=1)

                    print(FILESAVE.NON3)

                    cat("PM ", APPROXIMATION, " $(", NPSEUDO, ")$", " & $", pr1, "$ $(", pr2, ")$  \\\\\n\n", sep="")

                
              }
          }
            cat("\n\n", sep="")
      }
  }
}
