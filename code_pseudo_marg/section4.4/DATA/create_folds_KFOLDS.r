DATASETS = c("breast_cancer_wisconsin", "glass", "ionosphere", "pima", "thyroid", "usps")
## DATASETS = c("usps")
## DATASETS = c("pima")

NFOLDS = 5

set.seed(1)

for(DATASET in DATASETS)
{
filex = paste("clean/", DATASET, "_X.txt", sep="")
filey = paste("clean/", DATASET, "_y.txt", sep="")

X = as.matrix(read.table(filex, header=F))
y = as.matrix(read.table(filey, header=F))

n = dim(X)[1]
d = dim(X)[2]

ind.scramble = sample(c(1:n), n)
ind.folds = as.integer(seq(from=1, to=6-1e-6, length.out=n))

for(FOLD in c(1:NFOLDS))
  {
    ind.test = ind.scramble[which(ind.folds == FOLD)]
    ind.train = ind.scramble[which(ind.folds != FOLD)]
    tmpX.train = X[ind.train,]
    tmpy.train = y[ind.train]

    tmpX.test = X[ind.test,]
    tmpy.test = y[ind.test]

    mxi = sdxi = rep(0, d)
    for(i in 1:d)
      {
        mxi[i] = mean(tmpX.train[,i])
        sdxi[i] = sd(tmpX.train[,i])

        if(sdxi[i] == 0) browser()
        
        tmpX.train[,i] = tmpX.train[,i] - mxi[i]
        tmpX.train[,i] = tmpX.train[,i] / sdxi[i]

        tmpX.test[,i] = tmpX.test[,i] - mxi[i]
        tmpX.test[,i] = tmpX.test[,i] / sdxi[i]
      }
        
    print(sdxi)

    filetowritex = paste("folds/", DATASET, "_X_train_FOLD_", FOLD, ".txt", sep="")
    filetowritey = paste("folds/", DATASET, "_y_train_FOLD_", FOLD, ".txt", sep="")
    write.table(tmpX.train, row.names=F, col.names=F, quote=F, file=filetowritex)
    write.table(tmpy.train, row.names=F, col.names=F, quote=F, file=filetowritey)

    filetowritex = paste("folds/", DATASET, "_X_test_FOLD_", FOLD, ".txt", sep="")
    filetowritey = paste("folds/", DATASET, "_y_test_FOLD_", FOLD, ".txt", sep="")
    write.table(tmpX.test, row.names=F, col.names=F, quote=F, file=filetowritex)
    write.table(tmpy.test, row.names=F, col.names=F, quote=F, file=filetowritey)
  }
}
