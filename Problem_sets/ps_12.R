set.seed(191178)

#### Method 1 -> AR using U(-1, 1)

fx_ar <- function() {
    c <- 3/2
    accept <- 0
    try <- 0

    while (accept == 0) {
        try <- try + 1
        prop <- runif(1, min = -1, max = 1)
        U <- runif(1)
        ratio <- (1 - prop^2) 
        if (U < ratio) {
            return (c(prop, try))
        }
    }
}

N <- 1e5
samples <- numeric(length = N)
n.try <- numeric(length = N)

for (i in 1:N) {
    foo <- fx_ar()
    samples[i] <- foo[1]
    n.try[i] <- foo[2]
}

mean(samples)
mean(n.try)
par(mfrow = c(1,2))
hist(samples)
plot(density(samples))

#### Method 2 -> Ratio of Uniforms

N <- 1e6
samples <- matrix(0, nrow = N, ncol = 2)
done <- 0

while (done < N) {
    prop.u <- runif(1, min = 0, max = sqrt(3)/2)
    prop.v <- runif(1, min = -sqrt(3)/4, max = sqrt(3)/4)

    if (2*(prop.u^4) <= 3*(prop.u^2 - prop.v^2)) {
        done <- done + 1
        samples[done, ] = c(prop.u, prop.v)
    }
}

samp.rou <- samples[,2]/samples[,1]

mean(samp.rou)
par(mfrow = c(1,2))
hist(samp.rou)
plot(density(samp.rou))
