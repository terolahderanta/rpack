
# 03.04.2020

# Compare different clusterings with the corrected Rand index

# In particular, examine:

# 1. Do the capacity constraints improve the clustering when they are close to the "ground truth"

# 2. Does the outlier penalty improve the clustering

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
#library(parallel) # For parallel computing. rpack uses function "detectCores".
library(fpc) # For the corrected Rand index

# External functions

# These are added by Markku. These external functions generate data from gamma distribution.
# Weights are determined either by 1) chance 2 ) based on the distance from the distribution mean: the greater
# the Euclidean distance between the point and the population mean is, the larger the
# weight is. "simulate_unif_grid" is a sub-function of the main function "simulate_gamma_mixture".


source("CodeCollection/simulate_gamma_mixture.R")
source("CodeCollection/simulate_unif_grid.R")


# Initialize seed number:
#seed = Sys.time()
#seed = as.integer(seed)
#seed = seed %% 100000
seed = 36352
set.seed(seed)


# rpack clustering

# Compare clusterings (continuous setting):

# 1. Squared Euclidean, Lower and Upper limit, PLUS penalty

# 2. Squared Euclidean, Lower and Upper limit, NO penalty

# 3. Euclidean, Lower and Upper limits, PLUS penalty

# 4. Euclidean, Lower and Upper limit, NO penalty

# 5. Squared Euclidean NO Lower and Upper limit, NO penalty

# 6. Euclidean NO Lower and Upper limit, NO penalty

# 7. Squared Euclidean NO Lower and Upper limit, PLUS penalty

# 8. Euclidean NO Lower and Upper limit, PLUS penalty

###############################################3

# Outgroup parameter lambda, smaller value --> more outliers

lambda1 = 0.05 # For Squared Euclidean.

lambda2 = 0.2 # For Standard Euclidean

# In this artificial example, algorithm seems to return quite robust results for all
# tuning parameter values unless the penalty parameter is extremely small which returns
# numerous outliers

n_steps = 50 # Number of algorithm iterations

SimRounds = 50 # Nmb of simulation rounds

# Initialize the matrix where corrected Rand indexes are stored:

Results = matrix(0, SimRounds, 8)

# Number of clusters
k = 10

# Cluster sizes
n = c(20, 40, 60, 80, 50, 80, 60, 40, 20, 50) # uneven cluster sizes

n_out = 20 # nmb of outliers

for(i in 1:SimRounds){
  
  # First, generate 500 points from mixture of 10 gamma distributions.
  
  test_dat = simulate_gamma_mixture(n, k, n_out=n_out, out_scale = 5, scale_between_range = c(0, 1), outgroup_alpha = 0.4, place_on_grid = T,
                                    overlap_scale = 0.5)
  test_dat = test_dat$Y
  
  # Set the capacity constraints to true cluster sizes excluding outliers:
  
  pr_mean = round(sum(test_dat$w)/k)
  
  # Max radius for prior
  
  pr_width = 100
  
  # Lower und upper limit for cluster size
  
  L = (pr_mean - 10*pr_width)
  
  U = (pr_mean + 10*pr_width)
  
  # 1. Squared Euclidean, Lower and Upper limit, PLUS penalty
  
  clust1 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    weights = test_dat$w,
    N = n_steps,
    k = k,
    range = c(L, U),
    lambda = lambda1,
    place_to_point = FALSE
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 5] <- clust1$clusters
  
  # 2. Squared Euclidean, Lower and Upper limit, NO penalty
  
  # Alternaring algorithm
  clust2 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    weights = test_dat$w,
    N = n_steps,
    k = k,
    range = c(L, U),
    place_to_point = FALSE
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 6] <- clust2$clusters
  
  # 3. Euclidean, Lower and Upper limits, PLUS penalty
  
  clust3 <- alt_alg(
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
  
  test_dat[ , 7] <- clust3$clusters
  
  
  # 4. Euclidean, Lower and Upper limit, NO penalty
  
  # Alternaring algorithm
  clust4 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    weights = test_dat$w,
    N = n_steps,
    k = k,
    range = c(L, U),
    place_to_point = FALSE,
    d = euc_dist # Input for distance metric is in a form of a function.
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 8] <- clust4$clusters
  
  # 5. Squared Euclidean NO Lower and Upper limit, NO penalty
  
  # Alternaring algorithm
  clust5 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    weights = test_dat$w,
    N = n_steps,
    k = k,
    place_to_point = FALSE
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 9] <- clust5$clusters
  
  # 6. Euclidean NO Lower and Upper limit, NO penalty
  
  # Alternaring algorithm
  clust6 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    weights = test_dat$w,
    N = n_steps,
    k = k,
    place_to_point = FALSE,
    d = euc_dist
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 10] <- clust6$clusters
  
  # 7. Squared Euclidean NO Lower and Upper limit, PLUS penalty
  
  # Alternaring algorithm
  clust7 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    lambda = lambda1,
    weights = test_dat$w,
    N = n_steps,
    k = k,
    place_to_point = FALSE
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 11] <- clust7$clusters
  
  # 8. Euclidean NO Lower and Upper limit, PLUS penalty
  
  # Alternaring algorithm
  clust8 <- alt_alg(
    coords = dplyr::select(test_dat, x, y),
    lambda = lambda2,
    weights = test_dat$w,
    N = n_steps,
    k = k,
    place_to_point = FALSE,
    d = euc_dist
  )
  
  # Save the results of the clustering.
  
  test_dat[ , 12] <- clust8$clusters
  
  ####################################################
  
  # Determine the corrected Rand index:
  
  test_dat$orig_group = as.numeric(test_dat$orig_group)
  
  # Change the label of outgroup samples:
  
  for(j in 5:ncol(test_dat)){
    
    test_dat[ , j] = ifelse(test_dat[ , j] == 99, 11, test_dat[ , j])
    
  }
  
  d = dist(test_dat[, 1:2]) # Not needed for the Rand index but for the "cluster.stats".
  
  Results[i, 1] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,5])$corrected.rand
  Results[i, 2] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,6])$corrected.rand
  Results[i, 3] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,7])$corrected.rand
  Results[i, 4] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,8])$corrected.rand
  Results[i, 5] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,9])$corrected.rand
  Results[i, 6] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,10])$corrected.rand
  Results[i, 7] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,11])$corrected.rand
  Results[i, 8] = cluster.stats(d, test_dat[,4], alt.clustering = test_dat[,12])$corrected.rand
  
}

## Save the classification results

# 1. Squared Euclidean, Lower and Upper limit, PLUS penalty
# 2. Squared Euclidean, Lower and Upper limit, NO penalty
# 3. Euclidean, Lower and Upper limits, PLUS penalty
# 4. Euclidean, Lower and Upper limit, NO penalty
# 5. Squared Euclidean NO Lower and Upper limit, NO penalty
# 6. Euclidean NO Lower and Upper limit, NO penalty
# 7. Squared Euclidean NO Lower and Upper limit, PLUS penalty
# 8. Euclidean NO Lower and Upper limit, PLUS penalty

colnames(Results) = c("SqrdEuc_LU_plus_penalty",
                      "SqrdEuc_LU_no_penalty",
                      "Euc_LU_plus_penalty",
                      "Euc_LU_no_penalty",
                      "SqrdEuc_no_LU_no_penalty",
                      "Euc_no_LU_no_penalty",
                      "SqrdEuc_no_LU_plus_penalty",
                      "Euc_no_LU_plus_penalty")

par(cex.axis=0.5)
boxplot(Results)



#####

write.table(Results, "CodeCollection/ResultsSimulationFun/RandIndexResults.txt")
