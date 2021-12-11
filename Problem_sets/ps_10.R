set.seed(191178)

N <- 1e4
lambda <- 5
U <- runif(N)
samples <- -log(1-U)/lambda 

mean(samples)
par(mfrow=c(1,2))
hist(samples)
curve(dexp(x, rate = 5), 0, 1.5, col = 'red')
