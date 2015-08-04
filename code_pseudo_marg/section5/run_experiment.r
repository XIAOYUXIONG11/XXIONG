## DATASET = "breast_cancer_wisconsin"
## DATASET = "glass"
## DATASET = "ionosphere"
## DATASET = "pima"
## DATASET = "thyroid"
## DATASET = "usps"

## FOLD = 1

## METHOD = "SVM"
## METHOD = "EP_ML"
## METHOD = "EP_MCMC"
## METHOD = "EP_MCMC_PSEUDO"

## NTRAINING = 10

COVARIANCE.MODE = "D"
COVARIANCE.MODE = "H"

## ****************************** DATA

DATA = new.env()

filex = paste("~/PAMIDATA/folds_", NTRAINING, "/", DATASET, "_X_train_FOLD_", FOLD, ".txt", sep="")
filey = paste("~/PAMIDATA/folds_", NTRAINING, "/", DATASET, "_y_train_FOLD_", FOLD, ".txt", sep="")
DATA$X = as.matrix(read.table(filex, header=F))
DATA$y = read.table(filey, header=F)[[1]]
DATA$n = dim(DATA$X)[1]
DATA$d = dim(DATA$X)[2]

filex = paste("~/PAMIDATA/folds_", NTRAINING, "/", DATASET, "_X_test_FOLD_", FOLD, ".txt", sep="")
filey = paste("~/PAMIDATA/folds_", NTRAINING, "/", DATASET, "_y_test_FOLD_", FOLD, ".txt", sep="")
DATA$X.TEST = as.matrix(read.table(filex, header=F))
DATA$y.TEST = read.table(filey, header=F)[[1]]
DATA$ntest = dim(DATA$X.TEST)[1]

## ****************************** METHOD

if(METHOD == "SVM")
  {
    setwd("svm/")
    source("svm.r")
    setwd("..")
  }

if(METHOD == "EP_ML")
  {
    setwd("ep_ml/")
    source("ep_ml.r")
    setwd("..")
  }

if(METHOD == "EP_MCMC")
  {
    setwd("ep_mcmc/")
    source("ep_mcmc.r")
    setwd("..")
  }

if(METHOD == "EP_MCMC_PSEUDO")
  {
    setwd("ep_mcmc_pseudo/")
    source("ep_mcmc_pseudo.r")
    setwd("..")
  }

## ****************************** SAVE RESULTS

fileresults = paste("RESULTS_", NTRAINING, "/RESULTS_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, ".txt", sep="")
write.table(predictions.test, row.names=F, col.names=F, quote=F, file=fileresults)

fileresults.se = paste("RESULTS_", NTRAINING, "/RESULTS_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, "_STANDARD_ERROR.txt", sep="")
if(METHOD == "EP_MCMC") write.table(predictions.test.se, row.names=F, col.names=F, quote=F, file=fileresults.se)
if(METHOD == "EP_MCMC_PSEUDO") write.table(predictions.test.se, row.names=F, col.names=F, quote=F, file=fileresults.se)
