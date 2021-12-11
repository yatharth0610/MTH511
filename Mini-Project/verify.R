rm(list = ls())
load("191178.Rdata")

set.seed(191178)
dat <- read.table("https://dvats.github.io/assets/data/191178.txt")
X <- as.matrix(dat)

## Testing on complete dataset
test.error <- pred.loss(X.new = X, model = model)
test.error

## Testing on random 100 samples from the dataset
X.new <- X[sample(1:dim(X)[1], 100), ]
test.error <- pred.loss(X.new, model)
test.error


