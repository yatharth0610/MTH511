set.seed(191178)

N <- 1e4
lambda <- 1
alpha <- 1.5
U <- runif(N)
samples <- (-(log(1-U)/lambda))^(1/alpha) 

mean(samples)
hist(samples)
