set.seed(191178)

N <- 1e6
a <- 4
samples <- matrix(0, nrow = N, ncol = 2)
done <- 0

max_u <- 1 / sqrt(1 - exp(-a))
max_v <- 2/(exp(1) * sqrt(1 - exp(-a)))

while(done < N) {
    prop.u <- runif(1, min = 0, max = max_u)
    prop.v <- runif(1, min = 0, max = max_v)

    if (prop.v < a*prop.u && (prop.u)^2 <= (exp(-prop.v/prop.u))/(1 - exp(-a))) {
        done <- done + 1
        samples[done, ] = c(prop.u, prop.v)
    }

}

samp.rou <- samples[,2]/samples[,1]

mean(samp.rou)
par(mfrow = c(1,2))
hist(samp.rou)
plot(density(samp.rou), col = 'red', xlim = c(0, 4))
