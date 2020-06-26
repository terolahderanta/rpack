
# 31.03.2020

# The intersection of location-allocation and partitional clustering

# Artificial data example

# Don't use parallel computing because Windows Defender Firewall will just block it...

# Libraries

#First, load the necessary libraries. In addition to the `rpack` package, we use 
#`tidyverse` and `dplyr` for sample data manipulation and  `ggplot2` for plotting.

library(rpack)
library(tidyverse)
library(LaplacesDemon)
library(Matrix)
library(plotly)
library(Gmedian)
library(lcmix) # Added by Markku. Package can be found at R-Forge
library(parallel) # For parallel computing. rpack uses function "detectCores".
c_col = c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
          "chocolate","darkred","yellow3","darkgreen","bisque4","magenta",
          "royalblue","tomato4","steelblue1",
          "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
          "coral3","beige","brown4","indianred1","lightgreen","orchid")


# External functions

# These are added by Markku. These external functions generate data from gamma distribution.
# Weights are determined either by 1) chance 2 ) based on the distance from the distribution mean: the greater
# the Euclidean distance between the point and the population mean is, the larger the
# weight is. "simulate_unif_grid" is a sub-function of the main function "simulate_gamma_mixture".


source("CodeCollection/simulate_gamma_mixture.R")
source("CodeCollection/simulate_unif_grid.R")


# Simulated data 

#Let's set up some data to be clustered. To ease cluster overlap, clusters are placed on a grid. 

#Outliers are sampled uniformly on the cluster grid.

# Some simulated clusters have heavy weight points at edge of clusters (50/50 chance).** 


# Initialize seed number:
#seed = Sys.time()
#seed = as.integer(seed)
#seed = seed %% 100000
seed = 10376
set.seed(seed)
k = 10 # ten clusters
n = c(20, 40, 60, 80, 50, 80, 60, 40, 20, 50) # uneven cluster sizes
n_out = 20 # nmb of outliers
# Generating 500 points from mixture of 10 gamma distributions.
test_dat = simulate_gamma_mixture(n, k, n_out=n_out, out_scale = 5, scale_between_range = c(0, 1), outgroup_alpha = 0.4, place_on_grid = T,
overlap_scale = 0.5)
true_mu <- test_dat$mu_true
test_dat <- test_dat$Y
# Ids in interactive plot
id <-  1:nrow(test_dat)
plot_sim <- ggplot(data = test_dat, aes(x = x, y = y, size = w, label = id)) +
geom_point() +
scale_size(range = c(2, 6)) +  # Scale objects sizes
guides(
color = guide_legend(        # Point size in legend
override.aes = list(size=5)
)
) +
labs(x = "x", y = "y", title = "Unclustered data") +
theme(
legend.position = "right",            # Legend position and removing ticks from axis
axis.text.x = ggplot2::element_blank(),
axis.text.y = ggplot2::element_blank(),
axis.ticks = ggplot2::element_blank()
)
ggplotly(plot_sim, tooltip = c("id", "w"))


# rpack clustering

#The following steps are repeated in the future:

#1. Cluster the data into $k=10$ clusters.

#2. Set uniform prior for the cluster weights. In other words, we let the total weights of clusters vary
#uniformly within a given range. We set the range as $sum(w)/k - 1000; sum(w)/k + 1000$.

#3. Add outgroup for the data points. Outgroup consists of points that can be seen as outliers from other points and are not allocated to any of the clusters.


## Clustering with uniform prior: Squared Euclidean distance

#Calculate the distance between points as **squared** Euclidean distance. Squared Euclidean distance forms
#more spherical-like clusters(?)

### Continuous setting

#While clustering the data and cluster heads are defined, **cluster heads are not chosen from the data points.**
#This constraint problem formulation is referred to as a **continuous** location-allocation problem.


# Number of clusters
k <- 10

# Mean
pr_mean <- round(sum(test_dat$w) / k)

# Max radius for prior
pr_width <- 100

# Lower und upper limit for cluster size
L <- (pr_mean - 10*pr_width)
U <- (pr_mean + 10*pr_width)

# Outgroup parameter lambda, smaller value --> more outliers
lambda1 <- 0.05

# Number of iterations
n_steps = 50

# Alternaring algorithm
clust1_continuous <- alt_alg(
  coords = dplyr::select(test_dat, x, y),
  weights = test_dat$w,
  N = n_steps,
  k = k,
  range = c(L, U),
  lambda = lambda1,
  place_to_point = FALSE
)

# Save the results of the clustering.

test_dat$SqrE_limits_continuous <- clust1_continuous$clusters
mu_sqrE_limits_continuous <- clust1_continuous$centers


# Plot the clusters with `plot_clusters`.


plot_clusters(
coords = test_dat[,1:2],
weights = test_dat$w,
clusters = clust1_continuous$clusters,
centers = clust1_continuous$centers,
title = paste("Capacitated clustering, k = ", k, " (squared Euclidean, continuous)", sep = ""),
subtitle = paste("Uniform prior in [", L, ", ", U, "] on cluster sizes. Heavy weights at some cluster edges", sep = "")
)

# Squared Euclidean, continuous, NO lower and upper limits:

# Alternaring algorithm
clust2_continuous <- alt_alg(
  coords = dplyr::select(test_dat, x, y),
  weights = test_dat$w,
  N = n_steps,
  k = k,
  lambda = lambda1,
  place_to_point = FALSE
)

# Save the results of the clustering.

test_dat$SqrE_nolimits_continuous <- clust2_continuous$clusters
mu_sqrE_nolimits_continuous <- clust2_continuous$centers


# Plot the clusters with `plot_clusters`.


plot_clusters(
  coords = test_dat[,1:2],
  weights = test_dat$w,
  clusters = clust2_continuous$clusters,
  centers = clust2_continuous$centers,
  title = paste("Capacitated clustering, k = ", k, " (squared Euclidean, continuous, no limits)", sep = "")
)


## Clustering with different distance metric: Euclidean distance

# Once again, cluster the data into $k=10$ clusters, but calculate the distance between points as **Euclidean distance instead of default squared Euclidean distance**. Euclidean distance measure minimizes the physical distance between points and the center in a cluster, while squared Euclidean distance forms more spherical-like clusters.

### Continuous setting, lower and upper limits, NO outliers:


clust3_continuous <- alt_alg(
  coords = dplyr::select(test_dat, x, y),
  weights = test_dat$w,
  N = n_steps,
  k = k,
  range = c(L, U),
  place_to_point = FALSE,
  d = euc_dist # Input for distance metric is in a form of a function.  
)

# Save the results of the clustering.

test_dat$E_no_penalty_continuous <- clust3_continuous$clusters
mu_E_no_penalty_continuous <- clust3_continuous$centers


# Plot the clusters.


plot_clusters(
coords = test_dat[,1:2],
weights = test_dat$w,
clusters = clust3_continuous$clusters,
centers = clust3_continuous$centers,
title = paste("Capacitated clustering, k = " ,k, " (Euclidean, continuous, no outlier penalty)", sep = ""),
subtitle = paste("Uniform prior in [", L, ", ", U, "] on cluster sizes. Heavy weights at some cluster edges", sep = "")
)


# Continuous setting, Euclidean, outlier penalty


# Outgroup parameter lambda, smaller value --> more outliers
lambda2 <- 0.2

# Alternaring algorithm
clust4_continuous <- alt_alg(
  coords = dplyr::select(test_dat, x, y),
  weights = test_dat$w,
  N = n_steps,
  k = k,
  lambda = lambda2,
  range = c(L, U),
  place_to_point = FALSE,
  d = euc_dist # Input for distance metric is in a form of a function.
)

# Save the results of the clustering.

test_dat$E_limits_continuous <- clust4_continuous$clusters
mu_E_limits_continuous <- clust4_continuous$centers


#Plot the clusters.

plot_clusters(
coords = test_dat[,1:2],
weights = test_dat$w,
clusters = clust4_continuous$clusters,
centers = clust4_continuous$centers,
title = paste("Euclidean, continuous with penalty, k = ", k, sep = ""),
subtitle = paste("Uniform prior in [", L, ", ", U, "] on cluster sizes. Heavy weights at some cluster edges", sep = "")
)

## Save the allocation results

Centers = cbind(mu_sqrE_limits_continuous, mu_sqrE_nolimits_continuous,
                mu_E_no_penalty_continuous, mu_E_limits_continuous)

colnames(Centers) = c("SqrdEuc_continuous_x", "SqrdEuc_continuous_y",
                      "SqrdEuc_continuous_no_limit_x", "SqrdEuc_continuous_no_limit_y",
                      "Euc_continuous_no_penalty_x", "Euc_continuous_no_penalty_y",
                      "Euc_continuous_x", "Euc_continuous_y")

write.table(test_dat, "CodeCollection/ResultsSimulationFun/test_data.txt")
write.table(Centers, "CodeCollection/ResultsSimulationFun/Cluster_centers.txt", col.names = T)


## A bonus: the ground truth plot

# Plot the 'ground truth' clusters with `plot_clusters`.

k <- length(table(test_dat$orig_group))

trueC = plot_clusters(
 coords = test_dat[,1:2],
 weights = test_dat$w,
 clusters = test_dat$orig_group,
 centers = true_mu,
 title = paste("Ground truth clusters, k = " ,k, sep = ""),
 subtitle = "Heavy weights at some cluster edges"
) + theme(axis.text.x = element_text(angle=0),
          axis.text.y = element_text(angle = 0)) +
  coord_fixed()

trueC
