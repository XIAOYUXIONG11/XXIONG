## Code to generate the results presented in tables 1 and 3
## Comparison of the effect of the LA and EP approximations and number of importance samples in the efficiency of the sampling of hyper-parameters in probit classification

## Remove comments according to the parameters you want to set

## APPROACH = c("AA", "SURR", "PM")

## SEED = 1
## APPROXIMATION = "EP"

## NPSEUDO = 64

## NTRAINING = 50
## DTRAINING = 2

## COVARIANCE.MODE = "D"
## COVARIANCE.MODE = "H"

N = NTRAINING
D = DTRAINING

source("../functions/PRIORS.r")
source("../functions/GP.FUNCTIONS.r")


## ****************************** DATA

source("../functions/GENERATE.DATA.r")

## ****************************** METHOD

if(APPROACH == "PM") {
  ADAPTIVE = F

  setwd("mcmc_pseudo/")
  source("mcmc_pseudo.r")
  setwd("..")
}

if(APPROACH == "SURR") {
  ADAPTIVE = T

  setwd("surr/")
  source("gp.sampling.r")
  setwd("..")
}

if(APPROACH == "AA") {
  ADAPTIVE = T

  setwd("aa/")
  source("gp.sampling.r")
  setwd("..")
}

