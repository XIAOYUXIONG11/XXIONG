for(PLOT.TYPE in c("TRACE", "AUTOCORR", "CONV")) {

if(PLOT.TYPE == "TRACE") {
  main.plot.type="Trace plot"
  pdf.options(width=7,height=12,pointsize=28)
  ps.options(width=7,height=12,paper="special",horizontal=F, pointsize=28)
}

if(PLOT.TYPE == "AUTOCORR") {
  main.plot.type="Autocorr."
  pdf.options(width=3.5,height=12,pointsize=28)
  ps.options(width=3.5,height=12,paper="special",horizontal=F, pointsize=28)
}

if(PLOT.TYPE == "CONV") {
  main.plot.type="PSRF"
  pdf.options(width=3.5,height=12,pointsize=28)
  ps.options(width=3.5,height=12,paper="special",horizontal=F, pointsize=28)
}

library(coda)

DIR.RESULTS = "./"

DATASET.VECT = c("breast_cancer_wisconsin", "pima", "abalone")
DATASET.NAMES = c("Breast", "Pima", "Abalone")

FOLD = 1

METHOD.VECT = c("MCMC_AA", "MCMC_SURR", "MCMC_PSEUDO")
NPSEUDO.VECT = c(1)
COVARIANCE.MODE.VECT = c("H")

NSEEDS = 10
NPSEUDO = 1

for(DATASET in DATASET.VECT) {

  for(COVARIANCE.MODE in COVARIANCE.MODE.VECT) {

    if(PLOT.TYPE == "TRACE") {
      fileplot.pdf = paste("PLOT_TRACE_", DATASET, "_", COVARIANCE.MODE, ".pdf", sep="")
      fileplot.eps = paste("PLOT_TRACE_", DATASET, "_", COVARIANCE.MODE, ".eps", sep="")
    }
    if(PLOT.TYPE == "AUTOCORR") {
      fileplot.pdf = paste("PLOT_AUTOCORR_", DATASET, "_", COVARIANCE.MODE, ".pdf", sep="")
      fileplot.eps = paste("PLOT_AUTOCORR_", DATASET, "_", COVARIANCE.MODE, ".eps", sep="")
    }
    if(PLOT.TYPE == "CONV") {
      fileplot.pdf = paste("PLOT_CONV_", DATASET, "_", COVARIANCE.MODE, ".pdf", sep="")
      fileplot.eps = paste("PLOT_CONV_", DATASET, "_", COVARIANCE.MODE, ".eps", sep="")
    }

##     pdf(fileplot.pdf)
    postscript(fileplot.eps)
    par(mfrow=c(3,1))
    
    for(METHOD in METHOD.VECT) {
      
      zz.acceptance = rep(0, NSEEDS)
      zz.ess = rep(0, NSEEDS)
      
      SEED = 1
      if(METHOD == "MCMC_PSEUDO") {
        FILESAVECHAINS.psi.tau = paste(DIR.RESULTS, "mcmc_pseudo/CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
        FILESAVECHAINS.psi.sigma = paste(DIR.RESULTS, "mcmc_pseudo/CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
      }
      if(METHOD == "MCMC_AA") {
        FILESAVECHAINS.psi.tau = paste(DIR.RESULTS, "mcmc_aa/CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
        FILESAVECHAINS.psi.sigma = paste(DIR.RESULTS, "mcmc_aa/CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
      }
      if(METHOD == "MCMC_SURR") {
        FILESAVECHAINS.psi.tau = paste(DIR.RESULTS, "mcmc_surr/CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
        FILESAVECHAINS.psi.sigma = paste(DIR.RESULTS, "mcmc_surr/CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
      }
      
      
      NLS = as.integer(system(paste("head -1", FILESAVECHAINS.psi.tau, "| wc -w"), intern=T))
      
      tabellone = array(0, dim=c(12000, NSEEDS, NLS+1))
      
      for(SEED in 1:NSEEDS) {
      
        if(METHOD == "MCMC_PSEUDO") {
          main.method = "PM"
          FILESAVECHAINS.psi.tau = paste(DIR.RESULTS, "mcmc_pseudo/CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
          FILESAVECHAINS.psi.sigma = paste(DIR.RESULTS, "mcmc_pseudo/CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_NPSEUDO_", NPSEUDO, "_COV_", COVARIANCE.MODE, "_SEED_", SEED, ".txt", sep="")
        }
        if(METHOD == "MCMC_AA") {
          main.method = "AA"
          FILESAVECHAINS.psi.tau = paste(DIR.RESULTS, "mcmc_aa/CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
          FILESAVECHAINS.psi.sigma = paste(DIR.RESULTS, "mcmc_aa/CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
        }
        if(METHOD == "MCMC_SURR") {
          main.method = "SURR"
          FILESAVECHAINS.psi.tau = paste(DIR.RESULTS, "mcmc_surr/CHAINS_PSITAU_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
          FILESAVECHAINS.psi.sigma = paste(DIR.RESULTS, "mcmc_surr/CHAINS_PSISIGMA_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_HYP_MH_SEED_", SEED, ".txt", sep="")
        }
        
        tabellone[,SEED,1:NLS] = as.matrix(read.table(FILESAVECHAINS.psi.tau))
        tabellone[,SEED,NLS+1] = as.matrix(read.table(FILESAVECHAINS.psi.sigma))
        
      }

      main = paste(main.method, main.plot.type, sep=" - ")
      
      if(PLOT.TYPE == "TRACE") {
        par("mar"=c(2.0,2.5,1.1,0.3), "las"=1, "mgp"=c(1.8,0.6,0))
	if(DATASET == "abalone") YLIM = c(-1.2, 3)
	if(DATASET == "breast_cancer_wisconsin") YLIM = c(1.2, 3)
	if(DATASET == "pima") YLIM = c(0.9, 2.5)
        plot(tabellone[-c(1:2000), 1, 1][c(1:1000)*10], type="l", ylab="", xlab="", lwd=2, main=main, ylim=YLIM)
      }
      if(PLOT.TYPE == "AUTOCORR") {
        par("mar"=c(2.0,2.5,1.1,0.3), "las"=1, "mgp"=c(1.8,0.6,0))
        tmp = acf(tabellone[-c(1:2000), 1, 1][c(1:1000)*10], plot=F)
        plot(tmp[[1]], type="h", lwd=2, xlab="", ylab="", main=main, ylim=c(-1.96/sqrt(1000),1))
        abline(h=0)
        abline(h=1.96/sqrt(1000), lty=2)
        abline(h=-1.96/sqrt(1000), lty=2)
      }      
      if(PLOT.TYPE == "CONV") {
        tmp = mcmc.list()
        for(i in 1:NSEEDS) tmp[[i]] = mcmc(tabellone[1:12000,i,1])
      
        par("mar"=c(2.0,2.5,1.1,0.3), "las"=1, "mgp"=c(1.8,0.6,0))
        zz = gelman.plot(tmp, autoburnin=F, auto.layout=F, ylim=c(1, 3), main=main) 
        plot(zz[[2]], zz[[1]][,,'median'], type="l", ylab="", xlab="", lwd=2, main=main, ylim=c(1,3))
        points(zz[[2]], zz[[1]][,,'97.5%'], type="l", lwd=2, lty=2, col=2)
	abline(h=1)
      }
      
    }
    dev.off()
  }
}



}
