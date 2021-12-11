set.seed(191178)

#############
#### (a) part
#############

tnorm_ar <- function(a) {
    accept <- 0
    try <- 0
    c <- 1/(pnorm(a, mean = 0, sd = 1) - pnorm(-a, mean = 0, sd = 1)) + 0.01

    while (accept == 0) {
        try <- try + 1
        U <- runif(1)
        prop <- rnorm(1, mean = 0, sd = 1)

        if (prop^2 <= a^2) {
            ratio <- 1/c 
            if (U < ratio) {
                accept <- 1
                return (c(prop, try))
            }
        } 
    }
}

N <- 1e4
samples <- numeric(length = N)
n.try <- numeric(length = N)

### Performing for a = 4
for (i in 1:N) {
    foo <- tnorm_ar(4)
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
hist(samples)
mean(n.try)

samples <- numeric(length = N)
n.try <- numeric(length = N)

### Performing for a = 1
for (i in 1:N) {
    foo <- tnorm_ar(1)
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
hist(samples)
mean(n.try)

#############
#### (b) part
#############

## For p = 3

p <- 3
N <- 1e4
samples <- matrix(0, nrow = N, ncol = p)
done <- 0

# For a = 4
a <- 4
while (done < N) {
    prop <- rnorm(p, 0, 1)
    accept <- 1
    for (i in 1:p) {
        if (prop[i]^2 > a^2) {
            accept <- 0
            break
        }
    }

    if (accept) {
        done <- done + 1
        samples[done, ] <- prop
    }
}

par(mfrow = c(2,2))
plot(samples[,1], samples[,2], asp = 1)
plot(samples[,2], samples[,3], asp = 1)
plot(samples[,3], samples[,1], asp = 1)
plot(density(samples[,1]), main = "Marginal density for X1")

p <- 3
N <- 1e5
samples <- matrix(0, nrow = N, ncol = p)
done <- 0

# For a = 1
a <- 1
while (done < N) {
    prop <- rnorm(p, 0, 1)
    accept <- 1
    for (i in 1:p) {
        if (prop[i]^2 > a^2) {
            accept <- 0
            break
        }
    }

    if (accept) {
        done <- done + 1
        samples[done, ] <- prop
    }
}

par(mfrow = c(2,2))
plot(samples[,1], samples[,2], asp = 1)
plot(samples[,2], samples[,3], asp = 1)
plot(samples[,3], samples[,1], asp = 1)
plot(density(samples[,1]), main = "Marginal density for X1")