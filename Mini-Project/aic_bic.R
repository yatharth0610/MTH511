###########################################
## EM algorithm to fit mixture of Gaussians
## to multivariate data (old faithful)
##
## Cross-validation and AIC model selection
###########################################
library(mvtnorm)
library(ellipse)
set.seed(10)

data = read.table("https://dvats.github.io/assets/data/191178.txt")

# calculate negative log-likelihood
# of mixture of multivariate normal
log_like = function(X, pi.list, mu.list, Sigma.list, C)
{
  foo = 0
  for(c in 1:C)
  {
    foo = foo + pi.list[[c]]*dmvnorm(X, mean = mu.list[[c]], sigma = Sigma.list[[c]])
  }
  return(-sum(log(foo)))
}

# X is the data
# C is the number of clusters
GLMMforC = function(X, C, tol = 1e-3, maxit = 1e3)
{
  n = dim(X)[1]
  p = dim(X)[2]
  ######## Starting values ###################
  ## pi are equally split over C
  pi.list = rep(1/C, C)
	pdError = 1
	while(pdError){
		mu = list()  # list of all the means
		Sigma = list()  # list of all variances

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
			mu[[i]] =  rnorm(p, sd = 3) + colMeans(X)

			# another way to choose the starting means
			# mu[[i]] =  rnorm(p, sd = apply(X, 2, sd)) + colMeans(X)

			Sigma[[i]] = var(X)
		}
		# Choosing good starting values is important since
		# The GMM likelihood is not concave, so the algorithm
		# may converge to a local optima.

		######## EM algorithm steps ###################

		iter = 0
		diff = 100
		old.mu = mu
		old.Sigma = Sigma
		old.pi = pi.list
		epsilon = 1e-05

		Ep = matrix(0, nrow = n, ncol = C)  # gamma_{i,c}
		save.loglike = 0
		while((diff > tol) && (iter < maxit) )
		{
			iter = iter + 1
			## E step: find gammas
			for(c in 1:C)
			{
				Ep[ ,c] = pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
			}
			Ep = Ep/rowSums(Ep)

			### M-step
			pi.list = colMeans(Ep)
			for(c in 1:C)
			{
				mu[[c]] = colSums(Ep[ ,c] * X )/sum(Ep[,c])
			}

			for(c in 1:C)
			{
				foo = 0
				for(i in 1:n)
				{
					foo = foo + (X[i, ] - mu[[c]]) %*% t(X[i, ] - mu[[c]]) * Ep[i,c]
				}
				Sigma[[c]] = foo/sum(Ep[,c])

				if(sum(Ep[, c]) == 0 || min(eigen(Sigma[[c]])$values) <=0)
				{
					# Below is to ensure the estimator is positive definite
					# otherwise next iteration gamma_i,c,k cannot be calculated
					Sigma[[c]] = Sigma[[c]] + diag(epsilon, p)
					print("Matrix not positive-definite")
				}
			}
			save.loglike = c(save.loglike, log_like(X = X, pi.list = pi.list, mu.list = mu, Sigma.list = Sigma,  C = C))
			# Difference in the log-likelihoods as the difference criterion
			diff = abs(save.loglike[iter+1] - save.loglike[iter])

			old.mu = mu
			old.Sigma = Sigma
			old.pi = pi.list
		}
		pdError = 0

		# Final allocation updates
		for(c in 1:C)
		{
			Ep[ ,c] = pi.list[c]*apply(X, 1, dmvnorm , mu[[c]], Sigma[[c]])
		}
		Ep = Ep/rowSums(Ep)
	}
  return(list("pi" = pi.list, "mu" = mu, "Sigma" = Sigma, "Ep" = Ep, "log.like" = tail(save.loglike,1)))
}



# My data set
X = as.matrix(data)
n = dim(X)[1]
dim(X)
# Testing whether 2-7 classes are needed
potC = 3:6

######################################
##### Model selection via AIC #######
######################################

# aic = function(X, pi.list, mu.list, Sigma.list, C)
# {
#   nlike = log_like(X, pi.list, mu.list, Sigma.list, C)
#   rtn = 2*nlike + 2* (15*C-1) # No. of params  = 4C + C-1 + 10C
#   return(rtn)
# }

bic = function(X, pi.list, mu.list, Sigma.list, C)
{
  n = dim(X)[1]
  nlike = log_like(X, pi.list, mu.list, Sigma.list, C)
  rtn = 2*nlike + log(n)* (15*C-1) # No. of params  = 3*C -  1
  return(rtn)
}

# aicLike = numeric(length = length(potC))
bicLike = numeric(length = length(potC))
reps = 6
model.save = list()
model = list()
for(c in 1:length(potC))
{
  print(c)
  for(r in 1:reps)
  {
    model.save[[r]] = GLMMforC(X = X, C = potC[c])
  }

  chosen.run = which.min(sapply(model.save, function(t) t$log.like))
  model[[c]] = model.save[[chosen.run]]
  # aicLike[c] = aic(X = X, pi.list = model[[c]]$pi, mu.list =model[[c]]$mu, Sigma.list = model[[c]]$Sigma, C = potC[c])
  bicLike[c] = bic(X = X, pi.list = model[[c]]$pi, mu.list = model[[c]]$mu, Sigma.list = model[[c]]$Sigma, C = potC[c])
}

# aicLike
bicLike

par(mfrow = c(1,1))
chosen = which.min(bicLike)
# allot = apply(model[[chosen]]$Ep, 1, which.max)  ## Final allotment of classification
# plot(X[,1], X[,2], col = allot, pch = 16, main = "AIC: C = ") # plot allotment

# ell = list()
# for(c in 1:potC[[chosen]])
# {
#   ell[[c]] = ellipse(model[[chosen]]$Sigma[[c]], centre = as.numeric(model[[chosen]]$mu[[c]]))
#   lines(ell[[c]], col = c)
# }

# chosen = 1
allot = apply(model[[chosen]]$Ep, 1, which.max)  ## Final allotment of classification
plot(X[,1], X[,2], col = allot, pch = 16, main = "BIC: C = 4 maybe") # plot allotment

ell = list()
for(c in 1:potC[[chosen]])
{
  ell[[c]] = ellipse(model[[chosen]]$Sigma[[c]], centre = as.numeric(model[[chosen]]$mu[[c]]))
  lines(ell[[c]], col = c)
}

bicLike
