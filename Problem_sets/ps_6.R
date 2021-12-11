### f(x; n, p) = (x + n - 1)C(x)p^n(1-p)^x
set.seed(191178)

##### 
## Using Inverse Transform
#####

nbinom_inv <- function(r, p) {
    accept <- 0
    U <- runif(1)
    curr <- 0
    sum <- 0

    while (accept == 0) {
        sum <- sum + dnbinom(curr, size = r, prob = p)
        if (U < sum) {
            accept <- 1
            return (curr)
        }
        curr <- curr + 1
    }
}

r <- 10
p <- 0.30

N <- 1e5
samples <- numeric(length = N)

for (i in 1:N) {
    samples[i] = nbinom_inv(r, p)
}

mean(samples)
hist(samples)

###########
### Using Accept-Reject (Geometric Proposal)
###########

fx.mean <- r * (1 - p) / p
p.star <- 1 / (1 + fx.mean)

x <- 0:2e3
fx <- dnbinom(x, size = 10, prob = 0.30)
qx <- dgeom(x, prob = p.star)
max(fx/qx)

nbinom_ar <- function(r, p, p.star, c) {
    accept <- 0
    try <- 0

    while (accept == 0) {
        try <- try + 1
        U <- runif(1)
        prop <- rgeom(1, p.star)

        ratio <- dnbinom(prop, r, p) / (c * dgeom(prop, p.star))

        if (U < ratio) {
            accept <- 1
            return (c(prop, try))
        }

    }
}

N <- 1e4
c <- 8
samples <- numeric(length = N)
n.try <- numeric(length = N)

for (i in 1:N) {
    foo <- nbinom_ar(r, p, p.star, c)
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
mean(n.try)
hist(samples)
