DATASETS = c("breast_cancer_wisconsin", "glass", "ionosphere", "pima", "thyroid", "usps")
## DATASETS = c("usps")
DATASETS = c("banknote")

DATASETS = "abalone"

set.seed(122341)

for(DATASET in DATASETS)
{
filex = paste("clean/", DATASET, "_X.txt", sep="")
filey = paste("clean/", DATASET, "_y.txt", sep="")

X = as.matrix(read.table(filex, header=F))
y = as.matrix(read.table(filey, header=F))

n = dim(X)[1]
d = dim(X)[2]

ind.scramble = sample(c(1:n), n)

ind.train = ind.scramble
tmpX.train = X[ind.train,]
tmpy.train = y[ind.train]


mxi = sdxi = rep(0, d)
for(i in 1:d) {
  mxi[i] = mean(tmpX.train[,i])
  sdxi[i] = sd(tmpX.train[,i])
  
  if(sdxi[i] == 0) browser()
  
  tmpX.train[,i] = tmpX.train[,i] - mxi[i]
  tmpX.train[,i] = tmpX.train[,i] / sdxi[i]
}
        
print(sdxi)

filetowritex = paste("folds_notest/", DATASET, "_X_train.txt", sep="")
filetowritey = paste("folds_notest/", DATASET, "_y_train.txt", sep="")
write.table(tmpX.train, row.names=F, col.names=F, quote=F, file=filetowritex)
write.table(tmpy.train, row.names=F, col.names=F, quote=F, file=filetowritey)
}
