############################
########## ps_1 ############
############################

set.seed(191178)
data(cars)
cars

X <- cbind(1, cars[,1])
y <- as.matrix(cars[,-1], ncol = 1)

w <- qr.solve(t(X) %*% X) %*% t(X) %*% y 
w_real <- lm(cars[,2] ~ cars[,1], data = cars)

############################
########## ps_2 ############
############################
fuel2001 <- read.csv("https://dvats.github.io/assets/fuel2001.csv",
row.names = 1)

X <- as.matrix(cbind(1, fuel2001[,-c(2)]), nrow = dim(X)[1], ncol = dim(X)[2])
y <- as.matrix(fuel2001[,2], ncol = 1)

w <- qr.solve(t(X) %*% X, tol = 1e-10) %*% t(X) %*% y 
w_real <- lm(FuelC ~ ., data = fuel2001)

############################
########## ps_4 ############
############################
set.seed(1)
n <- 50
p <- 5
sigma2.star <- 1/2
beta.star <- rnorm(p)

X <- cbind(1, matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1)))
y <- X %*% beta.star + rnorm(n, mean = 0, sd = sqrt(sigma2.star))

beta <- qr.solve(t(X) %*% X) %*% t(X) %*% y 
sigma2 <- (t(y - X %*% beta) %*% (y - X %*% beta)) / n 

beta_ridge <- qr.solve(t(X) %*% X + diag(0.01, nrow = p)) %*% t(X) %*% y 
beta_ridge

beta_ridge <- qr.solve(t(X) %*% X + diag(1, nrow = p)) %*% t(X) %*% y 
beta_ridge

beta_ridge <- qr.solve(t(X) %*% X + diag(10, nrow = p)) %*% t(X) %*% y 
beta_ridge

beta_ridge <- qr.solve(t(X) %*% X + diag(100, nrow = p)) %*% t(X) %*% y 
beta_ridge

############################
########## ps_5 ############
############################
y.bar <- sum(y)/n
c.y <- y - y.bar

X.minus.one <- X[, -1]
X.bar <- X.minus.one - matrix(1, nrow = n, ncol = 1) %*% ((matrix(1, nrow = 1, ncol = n) %*% X.minus.one) / n) 

beta.minus.one <- qr.solve(t(X.bar) %*% X.bar) %*% t(X.bar) %*% c.y
beta.one <- y.bar - ((matrix(1, nrow = 1, ncol = n) %*% X.minus.one) / n) %*% as.matrix(beta.minus.one, ncol = 1)

beta <- t(cbind(beta.one, t(beta.minus.one)))

############################
########## ps_8 ############
############################
library(pracma)

alpha <- 4 # true value of alpha
n <- 10 # actual data size 

# X_1, X_2, \dots, X_n
dat <- rgamma(n, shape = alpha, rate = 1) 

alpha_newton <- numeric()
epsilon <- 1e-8  #some tolerance level preset
alpha_newton[1] <- 4  #alpha_0
count <- 1
tol <- 100 # large number

while(tol > epsilon)
{
  count <- count + 1
  
  #first derivative
  f.prime <- -n*psi(k = 0, alpha_newton[count - 1]) + sum(log(dat))
  
  #second derivative
  f.dprime <- -n*psi(k = 1, alpha_newton[count - 1])
  alpha_newton[count] <- alpha_newton[count - 1] - f.prime/f.dprime
  tol <- abs(alpha_newton[count] - alpha_newton[count-1])
}
alpha_newton
alpha.grid <- seq(0, 10, length = 100)
log.like <- numeric(length = 100)
for(i in 1:100)
{
  log.like[i] <- sum(dgamma(dat, shape = alpha.grid[i], log = TRUE))
}
plot(alpha.grid, log.like, type = 'l', xlab = expression(alpha), ylab = "Log Likelihood")
abline(v = alpha, col = "red", lty = 2)
for(t in 1:count)
{
  points(alpha_newton[t], sum(dgamma(dat, shape = alpha_newton[t], log = TRUE)), pch = 16)
}
abline(v = tail(alpha_newton[count]), col = "blue", lty = 2)

############################
########## ps_11 ###########
############################
mu.star <- 5
sigma2.star <- 2

n <- 100
dat <- rnorm(n, mean = mu.star, sd = sqrt(sigma2.star))

## beta using NR
mu_newton <- numeric()
epsilon <- 1e-8
mu_newton[1] <- 50
count <- 1
tol <- 100

while(tol > epsilon) {
    count <- count + 1
    mu_newton[count] <- mu_newton[count - 1] + sum(dat - mu_newton[count - 1]) / n
    tol <- abs(mu_newton[count] - mu_newton[count - 1])
}

mu_newton
mu <- mu_newton[count]

## sigma using NR
sigma_newton <- numeric()
epsilon <- 1e-8
sigma_newton[1] <- 0.2
count <- 1
tol <- 100

while(tol > epsilon) {
    count <- count + 1
    foo <- sum((dat - mu)^2)
    sigma <- sigma_newton[count - 1]

    f.prime <- (-n / sigma) + (foo / (sigma^3))
    f.dprime <- (n / sigma^2) - ((3 * foo) / (sigma^4))

    sigma_newton[count] <- sigma - f.prime/f.dprime
    tol <- abs(sigma_newton[count] - sigma)
    print(tol)
}
sigma.grid <- seq(0.5, 10, length = 100)
log.like <- numeric(length = 100)

for(i in 1:100)
{
  log.like[i] <- -n*log(sigma.grid[i]) - (n/2) * (log(2 * pi)) - (1/(2*(sigma.grid[i])^2)) * sum((dat - mu.star)^2)
  print(log.like[i])
}

plot(sigma.grid, log.like, type = 'l', xlab = expression(alpha), ylab = "Log Likelihood")
abline(v = sqrt(sigma2.star), col = "red", lty = 2)
for(t in 1:count)
{
  points(sigma_newton[t], -n*log(sigma_newton[t]) - (n/2) * (log(2 * pi)) - (1/(2*(sigma_newton[t])^2)) * sum((dat - mu.star)^2), pch = 16)
}
abline(v = tail(alpha_newton[count]), col = "blue", lty = 2)
sigma_newton

############################
########## ps_14 ###########
############################

alpha.star <- 5
beta.star <- 4

n <- 100
dat <- rbeta(n, alpha.star, beta.star)

hessian <- function(X, coeff) {
    alpha <- coeff[1]
    beta <- coeff[2]
    mat <- matrix(0, nrow = 2, ncol = 2)
    mat[1, 1] <- psi(k = 1, alpha) - psi(k = 1, alpha + beta)
    mat[1, 2] <- -psi(k = 1, alpha + beta)
    mat[2, 1] <- mat[1, 2]
    mat[2, 2] <- psi(k = 1, beta) - psi(k = 1, alpha + beta)
    return (mat)
}

grad <- function(X, coeff) {
    alpha <- coeff[1]
    beta <- coeff[2]
    vec <- matrix(0, nrow = 2, ncol = 1)
    vec[1] <- vec[1] + psi(k = 0, alpha) - psi(k = 0, alpha + beta) - (sum(log(X)) / n)
    vec[2] <- vec[2] + psi(k = 0, beta) - psi(k = 0, alpha + beta) - (sum(log(1 - X)) / n)
    return (vec)
}

epsilon <- 1e-8
tol <- 100
theta <- matrix(6, ncol = 1, nrow = 2)
iter <- 0

while (tol > epsilon) {
    iter <- iter + 1
    if (iter %% 1 == 0) print(tol)
    new <- theta - solve(hessian(X = dat, coeff = theta)) %*% grad(X = dat, coeff = theta)
    tol <- sum((theta - new)^2)
    theta <- new 
}

theta

############################
########## ps_19 ###########
############################
y1 <- 125
y2 <- 18
y3 <- 20
y4 <- 34

tol <- 100
epsilon <- 1e-8
theta <- 0
iter <- 0

while (tol > epsilon) {
  iter <- iter + 1
  if (iter %% 2 == 1) print(tol)
  new <- (y4 + (y1 * theta) / (theta + 2)) / (y2 + y3 + y4 + (y1 * theta) / (theta + 2))
  tol <- abs(new - theta)
  theta <- new
}

theta 

############################
########## ps_20 ###########
############################
y1 <- 125
y2 <- 18
y3 <- 20
y4 <- 34

tol <- 100
epsilon <- 1e-8
theta <- 0
iter <- 0

while (tol > epsilon) {
  iter <- iter + 1
  if (iter %% 2 == 1) print(tol)
  K <- 100
  samp <- rbinom(y1, K, theta/(theta + 2))
  exp <- sum(samp) / K
  new <- (y4 + exp)/(y2 + y3 + y4 + exp)
  tol <- abs(new - theta)
  theta <- new 
}

theta
