set.seed(191178)

bern_f <- function(p) {
    draw <- rbinom(5, 1, p)
    prod <- draw[1]*draw[2]*draw[3]*(1 - draw[4])*(1 - draw[5])
    return (prod)
}

N <- 1e4
p <- 0.4
samples <- numeric(length = N)

for (i in 1:N) {
    samples[i] <- bern_f(p)
}

mean(samples)
