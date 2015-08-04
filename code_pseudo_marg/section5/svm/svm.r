## Probabilistic SVM - grid search the parameters that maximize the internal cross validation error on the training set 

library(e1071)

C.vect = 10^(seq(from=-4, to=4, length.out=15))
gamma.vect = 10^(seq(from=-4, to=4, length.out=15))

model = list(accuracies = -Inf)
for(hhh in 1:length(gamma.vect))
  {
    for(jjj in 1:length(C.vect))
      {
        tmp.model = svm(DATA$X, as.factor(DATA$y), kernel = "radial", cost = C.vect[jjj], gamma = gamma.vect[hhh], cross=5, probability=T)
        print(mean(tmp.model$accuracies))
        if(mean(tmp.model$accuracies) > mean(model$accuracies)) model = tmp.model
      }
  }
cat("\n")

tmp = predict(model, DATA$X.TEST, probability=T)
predictions.test = attributes(tmp)$probabilities[,'1']

print(model$cost)
print(model$gamma)
