### f(x) = pi*f_1 + (1-pi)*(nCx)*p^x*(1-p)^(n-x)
set.seed(191178)

ZI_binom <- function(n, p, p.pi) {
    U <- runif(1)

    if (U < p.pi) {
        return (0)
    }    
    else {
        draw <- rbinom(1, n, p)
        return (draw)
    }
}

N <- 1e4
n <- 10
p <- 0.3
p.pi <- 0.5
samples <- numeric(length = N)

for (i in 1:N) {
    samples[i] <- ZI_binom(n, p, p.pi)
}

mean(samples)
hist(samples)
