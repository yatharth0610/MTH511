set.seed(191178)

discrete_inv_transform <- function(rate = p) {
    curr <- 0
    accept <- 0
    count <- 0

    draw <- runif(1, min = 0, max = 1)

    while(accept == 0) {
        curr <- curr + p*((1-p)^count)
        count <- count + 1

        if (draw < curr) {
            return (count)
        }
    }

}

N <- 1e4
samples <- numeric(length = N)
p <- 0.3

for (i in 1:N) {
    samples[i] <- discrete_inv_transform(p)
}

mean(samples)
hist(samples)
