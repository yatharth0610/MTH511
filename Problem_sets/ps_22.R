rm(list = ls())

set.seed(191178)

multinorm <- function(mu, sigma, N = 1e3) {
    decomp <- eigen(sigma)
    p <- length(mu)

    ## check for psd
    e_values <- decomp$values
    is_psd <- 1
    for (i in 1:p) {
        if (e_values[i] < 0) {
            is_psd <- 0
            break
        }
    }
    if (is_psd) {
        sigma.sq <- decomp$vectors %*% diag(e_values, p) %*% solve(decomp$vectors)
        samples <- matrix(0, nrow = N, ncol = p)

        for (i in 1:N) {
            prop <- rnorm(p, 0, 1)
            samples[i, ] <- mu + sigma.sq %*% prop
        }
        return(samples)
    }
}

mu <- c(5, 0, -2)
sigma <- matrix(c(1, .9, -.3, .9, 1, .1, -.3, .1, 1), nrow = 3, ncol = 3)
draws <- multinorm(mu, sigma)

par(mfrow = c(3,1))
plot(density(draws[,1]), main = "Marginal density of X1")
plot(density(draws[,2]), main = "Marginal density of X2")
plot(density(draws[,3]), main = "Marginal density of X3")

par(mfrow = c(3,1))
plot(draws[,1], draws[,2], main = "Correlation 0.9")
plot(draws[,2], draws[,3], main = "Correlation 0.1")
plot(draws[,1], draws[,3], main = "Correlation -0.3")
