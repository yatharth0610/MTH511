set.seed(191178)

gamma_ar <- function(alpha, beta, lambda, c) {
    accept <- 0
    try <- 0

    while (accept == 0) {
        try <- try + 1
        prop <- rexp(1, lambda)
        U <- runif(1)

        ratio <- dgamma(prop, alpha, beta) / (c * dexp(prop, lambda))

        if (U < ratio) {
            accept <- 1
            return (c(prop, try))
        }
    }
}

N <- 1e4
alpha <- 3/2
beta <- 1
lambda <- 2/3
c <- dgamma(3/2, alpha, beta) / dexp(2/3, lambda)

samples <- numeric(length = N) 
n.try <- numeric(length = N)

for (i in 1:N) {
    foo <- gamma_ar(alpha, beta, lambda, c)
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
par(mfrow = c(1, 2))
plot(density(samples))
curve(dgamma(x, alpha, beta), 0, 20, col = 'red')
mean(n.try)
