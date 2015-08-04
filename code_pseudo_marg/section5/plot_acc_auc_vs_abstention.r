## Code to generate the plots in figures 3, 4, and 5

library(Hmisc)
library(ROCR)

ps.options(width=10,height=10,paper="special",horizontal=F, pointsize=32)

capacity = function(pred, lab)
  {
    acc.abst = matrix(0, 0, 2)
    auc.abst = matrix(0, 0, 2)

    thrs.vect = seq(0, 0.5, by=0.01)
    for(thrs in thrs.vect)
      {
        iii = which( abs(pred - 0.5) >= thrs )
        if(length(iii) == 0) break
        if(length(levels(as.factor(y.TEST[iii])))==1) break
        tmp = prediction(p.predictions[iii], y.TEST[iii])
        auc.abst = rbind(auc.abst, c(1 - length(iii)/length(y.TEST), performance(tmp, "auc")@y.values[[1]]))

        tmp = table(c(pred[iii]>0.5, F, T), c(y.TEST[iii], -1, 1)) - diag(2)
        acc = sum(diag(tmp)) / sum(tmp)
        acc.abst = rbind(acc.abst, c(1 - length(iii)/length(y.TEST), acc))
      }

    nnn = dim(auc.abst)[1]

    if(DATASET == "glass") main = "Glass"
    if(DATASET == "ionosphere") main = "Ionosphere"
    if(DATASET == "pima") main = "Pima"
    if(DATASET == "thyroid") main = "Thyroid"
    if(DATASET == "usps") main = "USPS"
    
    if(METHOD == "SVM") {
      fileplot.AUC = paste("PLOT_VS_ABST_AUC_", DATASET, "_N_", NTRAINING, "_FOLD_", FOLD, ".eps", sep="")

      postscript(fileplot.AUC)
      par("mar"=c(2.9,2.9,1.1,0.4), "las"=1, "mgp"=c(1.8,0.6,0))
      
      plot(auc.abst[,1], auc.abst[,2], type='l', xlim=c(0,1), ylim=c(0,1), xlab="Abstention", ylab="AUC", pch=20, lwd=4, main=main)
    }

    if(METHOD != "SVM") points(auc.abst[,1], auc.abst[,2], type='l', xlim=c(0,1), ylim=c(0,1), col=INDMETHOD, pch=20, lwd=4)
    if(METHOD == "EP_MCMC_PSEUDO") {
      legend(0.4, 0.4, legend=c("SVM", "EP ML", "MCMC EP", "MCMC PM EP"), col=c(1:4), lwd=4)
      dev.off()
    }
    
    nnn = dim(auc.abst)[1]
    tmp1 = diff(acc.abst[,2]) / 2 + acc.abst[1:(nnn-1),2]
    tmp2 = diff(acc.abst[,1])
    capacity1 = sum(tmp1 * tmp2)

    tmp1 = diff(auc.abst[,2]) / 2 + auc.abst[1:(nnn-1),2]
    tmp2 = diff(auc.abst[,1])
    capacity2 = sum(tmp1 * tmp2)

##     ## Divide area by the maximum area (which is equal to the maximum value of abstention)
    capacity1 = capacity1 / auc.abst[nnn,1]
    capacity2 = capacity2 / auc.abst[nnn,1]

    list(capacity1, capacity2)
  }


NFOLDS = 40

DATASET = "breast_cancer_wisconsin"
DATASET = "glass"
DATASET = "ionosphere"
DATASET = "pima"
DATASET = "thyroid"
DATASET = "usps"

METHODS = c("SVM", "EP_MAP", "EP_MCMC", "EP_MCMC_PSEUDO")
NMETHODS = length(METHODS)

COVARIANCE.MODE = "H"


matrix.results.acc = matrix.results.auc = array(0, c(NFOLDS, NMETHODS, 4))

  iii.indtraining = 1
  
for(NTRAINING in c(10, 20, 50, 100))
{
  cat("N training", NTRAINING, "\n")
  for(FOLD in 1:NFOLDS)
    {
      cat(FOLD, "  \r")
      fileylabels = paste("~/PAMIDATA/folds_", NTRAINING, "/", DATASET, "_y_test_FOLD_", FOLD, ".txt", sep="")
      y.TEST = read.table(fileylabels, header=F)[[1]]
      
      for(INDMETHOD in 1:NMETHODS)
        {
          METHOD = METHODS[INDMETHOD]
          
          fileresults = paste("RESULTS_", NTRAINING, "/RESULTS_", DATASET, "_FOLD_", FOLD, "_METHOD_", METHOD, "_COV_", COVARIANCE.MODE, ".txt", sep="")
          p.predictions = read.table(fileresults, header=F)[[1]]
          
          tmp = capacity(p.predictions, y.TEST)
          
          matrix.results.acc[FOLD, INDMETHOD, iii.indtraining] = tmp[[1]]
          matrix.results.auc[FOLD, INDMETHOD, iii.indtraining] = tmp[[2]]
        }
    }
  
  iii.indtraining = iii.indtraining + 1
}

print(DATASET)
print(COVARIANCE.MODE)



################################################## Accuracy

if(DATASET == "glass") { limiti = c(0.5,1); limitisub = c(0.92, 1) }
if(DATASET == "ionosphere") limiti = c(0.3,1)
if(DATASET == "pima") limiti = c(0.48,0.9)
if(DATASET == "thyroid") { limiti = c(0.5,1); limitisub = c(0.96, 1) }
if(DATASET == "usps") { limiti = c(0.55,1); limitisub = c(0.93, 1) }

if(DATASET == "glass") legend.coords = c(2,0.9)
if(DATASET == "ionosphere") legend.coords = c(2,0.65)
if(DATASET == "pima") legend.coords = c(2,0.64)
if(DATASET == "thyroid") legend.coords = c(2,0.85)
if(DATASET == "usps") legend.coords = c(2,0.8)


if(DATASET == "glass") main = "Glass"
if(DATASET == "ionosphere") main = "Ionosphere"
if(DATASET == "pima") main = "Pima"
if(DATASET == "thyroid") main = "Thyroid"
if(DATASET == "usps") main = "USPS"

boxwex = 0.14
boxpasso = 0.16
quanti = 4

n.vect = c(10, 20, 50, 100)

res = list()
for(i in 1:4)
  {
    res[[i]] = data.frame(matrix.results.acc[,i,])
    names(res[[i]]) = n.vect
  }

fileplot = paste("PLOT_ACC_", DATASET, ".eps", sep="")


postscript(fileplot)
par("mar"=c(2.9,2.9,1.1,0.4), "las"=1, "mgp"=c(1.8,0.6,0))
boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab="n", ylab="Capacity accuracy", border=0, pch="", ylim=limiti, main=main)
boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, xlab="n", border=NULL, pch="", add=T, axes=F)
for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="") }

if(!(DATASET %in% c("pima", "ionosphere")))
  {
    if(DATASET == "thyroid") {
      splot = function() {
        boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
        axis(4)
        boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, tck=F)
        for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
      }

      subplot(splot(), size=c(5, 5), 1.2, 0.9, vadj=1, hadj=0)
    }

    if(DATASET == "glass") {
      splot = function() {
        boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
        axis(4)
        boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, tck=F)
        for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
      }

      subplot(splot(), size=c(5, 4), 1.2, 0.83, vadj=1, hadj=0)
    }

    if(DATASET == "usps") {
      splot = function() {
        boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
        axis(4)
        boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, tck=F)
        for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
      }
##       splot = function() {
##         quanti = quanti-1
##         boxplot(res[[2]][,2:4], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
##         axis(4)
##         boxplot(res[[2]][,2:4], at=c(1:quanti) - boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, boxfill=2, tck=F)
##         for(i in 3:4) { boxplot(res[[i]][,2:4], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
##       }

      subplot(splot(), size=c(4.7, 4.7), 1.4, 0.88, vadj=1, hadj=0)
    }

  }

if(DATASET %in% c("pima", "ionosphere")) legend(legend.coords[1], legend.coords[2], legend=c("SVM", "EP ML", "MCMC EP", "MCMC PM EP"), fill=c(0, 2:4))
dev.off()



################################################## AUC


if(DATASET == "glass") { limiti = c(0.05,1); limitisub = c(0.95, 1) }
if(DATASET == "ionosphere") limiti = c(0.3,1)
if(DATASET == "pima") limiti = c(0.31,0.92)
if(DATASET == "thyroid") {limiti = c(0,1); limitisub = c(0.95, 1) }
if(DATASET == "usps") { limiti = c(0.25,1); limitisub = c(0.97, 1) }

if(DATASET == "glass") legend.coords = c(2,0.9)
if(DATASET == "ionosphere") legend.coords = c(2,0.65)
if(DATASET == "pima") legend.coords = c(2,0.55)
if(DATASET == "thyroid") legend.coords = c(2,0.7)
if(DATASET == "usps") legend.coords = c(2,0.8)

if(DATASET == "glass") main = "Glass"
if(DATASET == "ionosphere") main = "Ionosphere"
if(DATASET == "pima") main = "Pima"
if(DATASET == "thyroid") main = "Thyroid"
if(DATASET == "usps") main = "USPS"


boxwex = 0.14
boxpasso = 0.16
quanti = 4

n.vect = c(10, 20, 50, 100)

res = list()
for(i in 1:4)
  {
    res[[i]] = data.frame(matrix.results.auc[,,i])
    names(res[[i]]) = n.vect
  }

fileplot = paste("PLOT_AUC_", DATASET, ".eps", sep="")


postscript(fileplot)
par("mar"=c(2.9,2.9,1.1,0.4), "las"=1, "mgp"=c(1.8,0.6,0))
boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab="n", ylab="Capacity AUC", border=0, pch="", ylim=limiti, main=main)
boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, xlab="n", border=NULL, pch="", add=T, axes=F)
for(i in 2:4)
  {
    boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="")
  }

if(!(DATASET %in% c("pima", "ionosphere")))
  {
    if(DATASET == "thyroid")
      {
        splot = function() {
          boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
          axis(4)
          boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, tck=F)
          for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
        }
        
        subplot(splot(), size=c(5, 4.5), 1.2, 0.73, vadj=1, hadj=0)
      }

    if(DATASET == "glass")
      {
        splot = function() {
          boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
          axis(4)
          boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, tck=F)
          for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
        }
        
        subplot(splot(), size=c(5, 5), 1.2, 0.8, vadj=1, hadj=0)
      }

    if(DATASET == "usps")
      {
        splot = function() {
          boxplot(res[[1]], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
          axis(4)
          boxplot(res[[1]], at=c(1:quanti) - 2 * boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, tck=F)
          for(i in 2:4) { boxplot(res[[i]], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
        }
##         splot = function() {
##           quanti = quanti-1
##           boxplot(res[[2]][,2:4], at=c(1:quanti), boxwex=boxwex, xlab='', ylab='', border=0, pch="", ylim=limitisub, yaxt="n")
##           axis(4)
##           boxplot(res[[2]][,2:4], at=c(1:quanti) - boxpasso, boxwex=boxwex, border=NULL, pch="", add=T, axes=F, boxfill=2, tck=F)
##           for(i in 3:4) { boxplot(res[[i]][,2:4], at=c(1:quanti) + (-3+i) * boxpasso, add=T, fill=T, names=NULL, axes=F, boxfill=i, boxwex=boxwex, pch="", tck=F) }
##         }
        
        subplot(splot(), size=c(5, 4.5), 1.2, 0.79, vadj=1, hadj=0)
      }
    
  }

if(DATASET %in% c("pima", "ionosphere")) legend(legend.coords[1], legend.coords[2], legend=c("SVM", "EP ML", "MCMC EP", "MCMC PM EP"), fill=c(0, 2:4))
dev.off()
