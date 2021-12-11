rm(list = ls())

set.seed(1)
n <- 50
nu <- 5
y <- rt(n, df = nu)

fx <- function(x) {
    ans <- exp(-(x^2) / 2)
    
    for (i in 1:50) {
        ans <- ans * (1 / (1 + (y[i] - x)/nu)^((nu + 1)/2))
    }

    return(ans)
}

solve <- function(N = 1e4) {
    prop <- rnorm(N)

    weights <- fx(prop) / dnorm(prop, 0, 1)
    theta <- mean(prop * weights) / mean(weights) 

    return(theta)
}

r <- 1e3
theta_hat <- numeric(length = r)
N <- 1e4

for (i in 1:r) {
    theta_hat[i] <- solve(N)
}

mean(theta_hat)
par(mfrow = c(1,2))
plot(theta_hat)
plot(1:r, theta_hat, type = "l", col = "green")
N*var(theta_hat)
