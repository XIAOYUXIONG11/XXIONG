## Code to produce the results to generate the plots in Figures 3 and 4

## Remove comments appropriately


## DATASET = "breast_cancer_wisconsin"
## DATASET = "pima"
## DATASET = "abalone"

FOLD = 1

## METHOD = c("MCMC_PSEUDO", "MCMC_AA", "MCMC_SURR")

## COVARIANCE.MODE = c("H")

## SEED = 1

## ****************************** DATA

DATA = new.env()

filex = paste("DATA/folds_notest/", DATASET, "_X_train.txt", sep="")
filey = paste("DATA/folds_notest/", DATASET, "_y_train.txt", sep="")
DATA$X = as.matrix(read.table(filex, header=F))
DATA$y = read.table(filey, header=F)[[1]]
DATA$n = dim(DATA$X)[1]
DATA$d = dim(DATA$X)[2]

## ****************************** METHOD

if(METHOD == "MCMC_PSEUDO")
  {
    setwd("mcmc_pseudo/")
    source("mcmc_pseudo.r")
    setwd("..")
  }

if(METHOD == "MCMC_AA")
  {
    setwd("mcmc_aa/")
    source("mcmc_aa.r")
    setwd("..")
  }

if(METHOD == "MCMC_SURR")
  {
    setwd("mcmc_surr/")
    source("mcmc_surr.r")
    setwd("..")
  }

## ****************************** 
