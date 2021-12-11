set.seed(191178)

N <- 1e5
samp <- rnorm(N, mean = 0, sd = 1)
draws <- samp * dt(samp, 1) / dnorm(samp, mean = 0, sd = 1)

plot(1:N, cumsum(draws)/(1:N), type = 'l', xlab = "N", ylab = "Running average")
