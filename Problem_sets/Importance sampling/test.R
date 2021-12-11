rm(list = ls())

set.seed(191178)

marginal <- function(N = 1e3, y, a)
{
  samp <- rnorm(N)
  mean(dt(y - samp, df = a))
}
grid <- 5e2
y.seq <- seq(-5, 5, length = grid)

f.y <- numeric(length = grid)
for(i in 1:length(y.seq))
{
  f.y[i] <- marginal(N = 1e3, y = y.seq[i], a = 3)
}
plot(y.seq, f.y, type= 'l', xlab = "Y", ylab = "Estimated density")

#######################

N <- 1e4
r <- 1e3
samples <- numeric(length = r)

for (i in 1:r) {
  props <- rnorm(N, 0, 4)
  fx <- props * dnorm(props, 0, 1) / dnorm(props, 0, 4)
  samples[i] <- mean(fx)
}

mean(samples)
N * var(samples)

########################
 
c <- 3

tgamma <- function(alpha, beta) {
  N <- 1e4
  r <- 1e2
  theta <- numeric(length = r)

  for (i in 1:r) {
    samples <- numeric(length = N)
    done <- 0

    while (done < N) {
      prop <- rgamma(1, alpha, beta)
      if (prop > c) {
        done <- done + 1
        samples[done] <- prop
      }
    }

    theta[i] <- mean(samples)
  }

  return(c(mean(theta), var(theta)))
}

foo <- matrix(0, nrow = 3, ncol = 2)

foo[1, ] <- tgamma(4,2)
foo[2, ] <- tgamma(8,2)
foo[3, ] <- tgamma(5,2)

foo
#########################

