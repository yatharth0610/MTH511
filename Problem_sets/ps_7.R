set.seed(191178)

fx <- function(x, lambda, m) {
    sum <- 0

    for (i in 0:m) {
        sum <- sum + dpois(i, lambda)
    }

    return (dpois(x, lambda)/sum)
}

tpois_ar <- function(lambda, m, c) {
    accept <- 0
    try <- 0

    while (accept == 0) {
        try <- try + 1
        U <- runif(1)
        prop <- sample(1:m, 1)
        ratio <- 30*fx(prop, lambda, m) / c 

        if (U < ratio) {
            accept <- 1
            return (c(prop, try))
        }
    }

}

m <- 30
lambda <- 20

all_x <- 0:m 
all_c <- 30 * fx(all_x, lambda, m)
c <- max(all_c) + 0.01
c

N <- 1e4
samples <- numeric(length = N)
n.try <- numeric(length = N)

for (i in 1:N) {
    foo <- tpois_ar(lambda, m, c)
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
hist(samples)
mean(n.try)
