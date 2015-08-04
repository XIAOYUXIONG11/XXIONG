## Code to create folds with n = NTRAINING

NTRAINING = 10

DATASETS = c("breast_cancer_wisconsin", "glass", "ionosphere", "pima", "thyroid", "usps")

NFOLDS = 40

for(DATASET in DATASETS)
{
filex = paste("clean/", DATASET, "_X.txt", sep="")
filey = paste("clean/", DATASET, "_y.txt", sep="")

X = as.matrix(read.table(filex, header=F))
y = as.matrix(read.table(filey, header=F))

n = dim(X)[1]
d = dim(X)[2]

for(FOLD in c(1:NFOLDS))
  {
    set.seed(FOLD)
    
    ind.ones = which(y == 1)
    ind.minus.ones = which(y == -1)

    ind.train = c(sample(ind.ones, NTRAINING/2), sample(ind.minus.ones, NTRAINING/2))
    ind.train = ind.train[sample(c(1:NTRAINING), NTRAINING)]

    ind.test = c(1:n)[-ind.train]
    tmpX.train = X[ind.train,]
    tmpy.train = y[ind.train]

    tmpX.test = X[ind.test,]
    tmpy.test = y[ind.test]

    mxi = sdxi = rep(0, d)
    for(i in 1:d)
      {
        mxi[i] = mean(tmpX.train[,i])
        sdxi[i] = sd(tmpX.train[,i])

        if(sdxi[i] != 0) 
          {
            tmpX.train[,i] = tmpX.train[,i] - mxi[i]
            tmpX.train[,i] = tmpX.train[,i] / sdxi[i]

            tmpX.test[,i] = tmpX.test[,i] - mxi[i]
            tmpX.test[,i] = tmpX.test[,i] / sdxi[i]
          }
      }
        
    filetowritex = paste("folds_", NTRAINING, "/", DATASET, "_X_train_FOLD_", FOLD, ".txt", sep="")
    filetowritey = paste("folds_", NTRAINING, "/", DATASET, "_y_train_FOLD_", FOLD, ".txt", sep="")
    write.table(tmpX.train, row.names=F, col.names=F, quote=F, file=filetowritex)
    write.table(tmpy.train, row.names=F, col.names=F, quote=F, file=filetowritey)

    filetowritex = paste("folds_", NTRAINING, "/", DATASET, "_X_test_FOLD_", FOLD, ".txt", sep="")
    filetowritey = paste("folds_", NTRAINING, "/", DATASET, "_y_test_FOLD_", FOLD, ".txt", sep="")
    write.table(tmpX.test, row.names=F, col.names=F, quote=F, file=filetowritex)
    write.table(tmpy.test, row.names=F, col.names=F, quote=F, file=filetowritey)
  }
}
