set.seed(191178)

binom_ar <- function(n, p, p.star, c) {
    accept <- 0
    count <- 0

    while (accept == 0) {
        count <- count + 1
        unif_draw <- runif(1)
        prop <- rgeom(1, p.star)
        
        if (prop <= n) {
            temp <- dbinom(prop, n, p) / (c * dgeom(prop, p.star))

            if (unif_draw < temp) {
                accept <- 1
                return (c(prop, count))
            }
        }
    }
}

N <- 1e6
p <- 0.75
n <- 20
p.star <- 1 / (n*p + 1)
x <- 0:n
all_c <- choose(n,x) * (1-p)^(n - x) *(1-p.star)^(-x) * p^(x) / p.star
(c <- max(all_c) + 0.01)
samples <- numeric(length = N)
n.try <- numeric(length = N)

for (i in 1:N) {
    foo <- binom_ar(n, p, p.star, c)
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
mean(n.try)
plot(density(samples))
