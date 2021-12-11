set.seed(10)
library(mvtnorm)
library(ellipse)

dat <- read.table("https://dvats.github.io/assets/data/191178.txt")
X <- as.matrix(dat)
n = dim(X)[1]
# calculate negative log-likelihood
# of mixture of multivariate normal
log_like <- function(X, pi.list, mu.list, Sigma.list, C)
{
  foo <- 0
  for(c in 1:C)
  {
    foo <- foo + pi.list[[c]]*dmvnorm(X, mean = mu.list[[c]], sigma = Sigma.list[[c]])
  }
  return(-sum(log(foo)))
}

GLMMforC <- function(X, C, tol = 1e-5, maxit = 1e3)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  ######## Starting values ###################
  ## pi are equally split over C
  pi.list <- rep(1/C, C)
  
  mu <- list()  # list of all the means
  Sigma <- list()  # list of all variances
  
  # The means for each C cannot be the same, 
  # since then the three distributions overlap
  # Hence adding random noise to colMeans(X)
  for(i in 1:C)
  {
    # cannot start with the same mean
    # for each component. Thus adding
    # random noise to the starting means
    # this also encourages convergence to
    # different solutions
    # mu[[i]] <-  rnorm(p, sd = 3) + colMeans(X)

    # another way to choose the starting means
    mu[[i]] <-  rnorm(p, sd = apply(X, 2, sd)) + colMeans(X)   
     
    Sigma[[i]] <- var(X)
  }
  # Choosing good starting values is important since
  # The GMM likelihood is not concave, so the algorithm
  # may converge to a local optima.
  
  ######## EM algorithm steps ###################
  
  iter <- 0
  diff <- 100
  old.mu <- mu
  old.Sigma <- Sigma
  old.pi <- pi.list
  epsilon <- 1e-05
  
  Ep <- matrix(0, nrow = n, ncol = C)  # gamma_{i,c}
  save.loglike <- 0
  while((diff > tol) && (iter < maxit) )
  {
    iter <- iter + 1
    flag <- 0 
    ## E step: find gammas
    for(c in 1:C)
    {
      Ep[ ,c] <- pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
    }
    Ep <- Ep/rowSums(Ep)
    
    ### M-step
    pi.list <- colMeans(Ep)
    for(c in 1:C)
    {
      mu[[c]] <- colSums(Ep[ ,c] * X )/sum(Ep[,c])
    }
    
    for(c in 1:C)
    {
      foo <- 0
      for(i in 1:n)
      {
        foo <- foo + (X[i, ] - mu[[c]]) %*% t(X[i, ] - mu[[c]]) * Ep[i,c] 
      }
      Sigma[[c]] <- foo/sum(Ep[,c])

      if(sum(Ep[, c]) == 0 || min(eigen(Sigma[[c]])$values) <=0)
      {
        # Below is to ensure the estimator is positive definite
        # otherwise next iteration gamma_i,c,k cannot be calculated
        Sigma[[c]] <- Sigma[[c]] + diag(epsilon, p)
        print("Matrix not positive-definite")
        flag <- 1
        pi.list <- rep(1/C, C)
  
        mu <- list()  # list of all the means
        Sigma <- list()  # list of all variances
        
        # The means for each C cannot be the same, 
        # since then the three distributions overlap
        # Hence adding random noise to colMeans(X)
        for(i in 1:C)
        {
          # cannot start with the same mean
          # for each component. Thus adding
          # random noise to the starting means
          # this also encourages convergence to
          # different solutions
          # mu[[i]] <-  rnorm(p, sd = 3) + colMeans(X)

          # another way to choose the starting means
          mu[[i]] <-  rnorm(p, sd = apply(X, 2, sd)) + colMeans(X)   
          
          Sigma[[i]] <- var(X)
        }

        iter <- 0
        diff <- 100
        old.mu <- mu
        old.Sigma <- Sigma
        old.pi <- pi.list
        epsilon <- 1e-05

        Ep <- matrix(0, nrow = n, ncol = C)  # gamma_{i,c}
        save.loglike <- 0

        break

      }
    }
    
    if (flag) {
      next
    }

    save.loglike <- c(save.loglike, log_like(X = X, pi.list = pi.list, mu.list = mu, Sigma.list = Sigma,  C = C))
    # Difference in the log-likelihoods as the difference criterion
    diff <- abs(save.loglike[iter+1] - save.loglike[iter])
    
    old.mu <- mu
    old.Sigma <- Sigma
    old.pi <- pi.list
  }

  print(iter)
  if (iter == maxit) {
    print("Not converged yet!")
  }
  
  # Final allocation updates
  for(c in 1:C)
  {
    Ep[ ,c] <- pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
  }
  Ep <- Ep/rowSums(Ep)
  
  return(list("Clusters" = C, "pi" = pi.list, "mu" = mu, "Sigma" = Sigma, "log.like" = tail(save.loglike,1)))
}

## Trying with Cross Validation

################################################
# 5-fold Cross-validation

# permutation <- sample(1:n, replace = FALSE)
# K <- 5
# # Uneven folds, but that is ok
# test.index <- split(permutation, rep(1:K, length = n, each = n/K))

# # Testing whether 2-4 classes are needed
# potC <- 6:7
# CV.errorLike <- numeric(length = length(potC))

# # will run EM multiple times for each training data
# # since EM convergence to local minima. Setting these
# # reps = 7
# reps <- 5
# model.save <- list()
# for(c in 1:length(potC))
# {
#   foo3 <- 0
#   for(k in 1:K)
#   {
#     print(c(c,k))
#     X.train <- X[-test.index[[k]], ]
#     X.test <- X[test.index[[k]], ]
    
#     for(r in 1:reps)
#     {
#       model.save[[r]] <- GLMMforC(X = X.train, C = potC[c])
#     }
#     # which ever run is the lowest negative log-like
#     chosen.run <- which.min(sapply(model.save, function(t) t$log.like))
#     model <- model.save[[chosen.run]] 
#     foo3 <- foo3 + log_like(X = X.test, pi.list = model$pi, mu.list =model$mu, Sigma.list = model$Sigma,  C = potC[c])
#   }
#   CV.errorLike[c] <- foo3/n
# }

# CV.errorLike

## Trying with AIC-BIC losses

aic <- function(X, pi.list, mu.list, Sigma.list, C)
{
  nlike <- log_like(X, pi.list, mu.list, Sigma.list, C)
  rtn <- 2*nlike + 2* (15*C - 1)
  return(rtn)
}

bic <- function(X, pi.list, mu.list, Sigma.list, C)
{
  n <- dim(X)[1]
  nlike <- log_like(X, pi.list, mu.list, Sigma.list, C)
  rtn <- 2*nlike + log(n)* (15*C - 1) 
  return(rtn)
}

potC <- 5

aicLike <- numeric(length = length(potC))
bicLike <- numeric(length = length(potC))
reps <- 6
model.save <- list()
model <- list()
for(c in 1:length(potC))
{
  print(c)
  for(r in 1:reps)
  {
    model.save[[r]] <- GLMMforC(X = X, C = potC[c])
  } 
  
  chosen.run <- which.min(sapply(model.save, function(t) t$log.like))
  model[[c]] <- model.save[[chosen.run]]
  aicLike[c] <- aic(X = X, pi.list = model[[c]]$pi, mu.list =model[[c]]$mu, Sigma.list = model[[c]]$Sigma, C = potC[c])
  bicLike[c] <- bic(X = X, pi.list = model[[c]]$pi, mu.list =model[[c]]$mu, Sigma.list = model[[c]]$Sigma, C = potC[c])
}

aicLike 
bicLike  

best_model <- which.min(bicLike)
model <- model[[best_model]]

pred.loss <- function(X.new, model)
{
  library(mvtnorm)
  tot <- 0
  pi.list <- model$pi 
  mu.list <- model$mu
  Sigma.list <- model$Sigma

  for(c in 1:model$Clusters)
  {
    tot <- tot + pi.list[[c]]*dmvnorm(X.new, mean = mu.list[[c]], sigma = Sigma.list[[c]])
  }
  loss <- -sum(log(tot))
  return(loss)
}

# par(mfrow = c(1,1))
# allot <- apply(model$Ep, 1, which.max)  ## Final allotment of classification
# plot(X[,1], X[, 2], col = allot, pch = 16,
#      main = paste("Log-Like = ", round(model$log.like,3))) # plot allotment

# ell <- list()
# for(c in 1:model$Clusters)
# {
#   ell[[c]] <- ellipse(model$Sigma[[c]], centre = as.numeric(model$mu[[c]]))
#   lines(ell[[c]], col = c)
# }  

pred.loss(X, model)

save(pred.loss, model, file = "191178.Rdata")
plot(dat)
