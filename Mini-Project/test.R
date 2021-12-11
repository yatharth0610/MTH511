set.seed(10)
library(mvtnorm)
library(mclust)

dat <- read.table("https://dvats.github.io/assets/data/191178.txt")
X <- as.matrix(dat)
BIC <- mclustBIC(X)
mod1 <- Mclust(X, x = BIC)
summary(mod1, parameters = TRUE)
X <- as.matrix(scale(dat))
BIC <- mclustBIC(X)
mod2 <- Mclust(X, x = BIC)
summary(mod2, parameters = TRUE)

# plot(mod1, what = "classification")
# opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 10, criterion = "BIC", 
                               
#                                dist_mode = "maha_dist", seed_mode = "random_subset",
                               
#                                km_iter = 10, em_iter = 10, var_floor = 1e-10, 
                               
#                                plot_data = T)

# new_dat <- tsne(dat, perplexity = 400)
# plot(new_dat)


# dat <- read.csv("https://dvats.github.io/assets/battingbowling.csv")
# X <- apply(as.matrix(dat)[,c(2,3)], 2, as.numeric)
# opt_gmm = Optimal_Clusters_GMM(X, max_clusters = 7, criterion = "BIC", 
                               
#                                dist_mode = "maha_dist", seed_mode = "random_subset",
                               
#                                km_iter = 10, em_iter = 10, var_floor = 1e-10, 
                               
#                                plot_data = T)

# data(faithful)
# dat <- as.matrix(faithful)
# plot(dat, pch = 16)

# opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 7, criterion = "BIC", 
                               
#                                dist_mode = "maha_dist", seed_mode = "random_subset",
                               
#                                km_iter = 10, em_iter = 10, var_floor = 1e-10, 
                               
#                                plot_data = T)

