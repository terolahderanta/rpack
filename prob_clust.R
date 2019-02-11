# Probabilistic clustering algorithm where constraints on cluster sizes are given as a prior distribution. 

# Packages ----------------------------------------------------------------

library(lpSolveAPI)
library(lpSolve)
library(mvtnorm)
c_col = c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple","chocolate","darkgrey","darkred","yellow3","darkgreen","wheat3","magenta","palegreen2","violetred","seagreen2","tomato4","steelblue1","royalblue","seagreen4","orangered","darkblue","khaki3","lavender","deeppink2","coral3","beige","brown4","indianred1","lightgreen","orchid")


# Functions ---------------------------------------------------------------


getmode <- function(v) {
  # Returns the mode of vector v
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

prob_clust <- function(data, weights, k, init_mu, prior_cl_sizes, prior_prob){
  # Alternatinh algorithm for maximizing the joint density.
  # Parameters
  #   data: matrix or data.frame containing the data.
  #   weights: vector of weights for each data point.
  #   init_mu: parameters (locations) that define the k distributions.
  #   k: number of clusters.
  #   prior_cl_sizes: all the possible values for cluster sizes.
  #   prior_prob: corresponding probabilities (sum to 1).  
  # Returns
  #   clusters: new cluster allocations for each object in data
  #   mu: new cluster center locations
  
  # Number of objects in data
  n <- length(data[,1])
  
  # Weights must be integers
  weights <- round(weights)
  
  # Equally weighted data
  data_ew <- apply(data, MARGIN = 2, FUN = rep, times = weights)
  
  # Original id:s for each data point in equally weighted data
  id_ew = rep(1:n, times = weights)
  
  # Initial clusters
  clusters <- rep(0,n)
  
  #Cluster centers
  mu <- init_mu
  
  max_sim <- 50
  for (iter in 1:max_sim) {
    # Old mu is saved to check for convergence
    old_mu <- mu
    
    # Clusters in equally weighted data (M-step)
    clusters_ew <- prob_clust_mstep(data_ew, mu, k, prior_cl_sizes, prior_prob)
    
    # Converting clusters_ew to original data
    for (i in 1:n) {
      clusters[i] <- getmode(clusters_ew[id_ew == i])
    }
    
    # Updating cluster centers (E-step)
    mu <- prob_clust_estep(data, weights, clusters, k)
    
    print(paste("Iteration:",iter))
    
    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }
  return(list(clusters,mu))
}

prob_clust_estep <- function(data, weights, clusters, k){
  # Updates the parameters (centers) for each cluster.
  # Parameters
  #   data: matrix or data.frame containing the data.
  #   weights: vector of weights for each data point.
  #   clusters: vector of cluster assignments for each data point.
  #   k: number of clusters.
  # Returns
  #   mu: new cluster centers.
  
  # Matrix for cluster centers
  mu <- matrix(0,nrow = k, ncol = length(data[1,]))
  
  # Update mu given 
  for (i in 1:k) {
    mu[i,] <- apply(data[clusters == i,] * weights[clusters == i], 2, FUN=sum) / sum(weights[clusters == i])  
  }
  
  return(mu)
}

prob_clust_mstep <- function(data_ew, mu, k, prior_cl_sizes, prior_prob){
  # Updates cluster allocations by maximizing the joint log-likelihood (M-step).
  # Parameters
  #   data_ew: matrix or data.frame containing the data, where each object is considered to be equally weighted.
  #   mu: parameters (locations) that define the k distributions.
  #   k: number of clusters.
  #   prior_cl_sizes: all the possible values for cluster sizes.
  #   prior_prob: corresponding probabilities (sum to 1).
  # Returns
  #   new_allocation: new cluster allocations for each object in data_ew
  
  # Number of objects in data_ew
  n <- length(data_ew[,1])
  
  # Log priors for possible cluster sizes
  prior_log_prob <- log(prior_prob)
  
  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(mu, MARGIN = 1, FUN = dmvnorm, x = data_ew, sigma = diag(2), log = TRUE)
  
  # Number of possible cluster sizes
  A <- length(prior_cl_sizes)
  
  # New linear program model object containing, where the number of decision variables is n * k + k * A
  lp1 <- make.lp(nrow = 0, ncol = n * k + k * A)
  
  # Maximizing the joint log-likelihood
  lp.control(lp1, sense = "max")
  
  # Binary integer problem
  set.type(lp1, 1:(n*k + k*A), "binary")
  
  # Objective function in the optimization problem
  set.objfn(lp1, c(C, t(matrix(1, ncol = k, nrow = A) * prior_log_prob)))
  
  # Constraint 1: each point is assigned to exactly one cluster
  for (i in 1:n) {
    add.constraint(lp1, c(as.numeric(rep(1:n == i, k)), rep(0, k * A)), "=", 1)
  }
  
  # Constraint 2: each cluster has exactly one size
  for (i in 1:k) {
    add.constraint(lp1, c(rep(0, k * n), as.numeric(rep(1:k == i, A))), "=", 1)
  }
  
  # Constraint 3: each cluster has S_a points
  temp <- c(t(matrix(1, ncol = n, nrow = k) * 1:k))
  temp2 <- c(t(-matrix(1, ncol = k, nrow = A) * prior_cl_sizes))
  for (i in 1:k) {
    add.constraint(lp1, c(as.numeric(temp == i, n), temp2 * as.numeric(rep(1:k == i, A))), "=", 0)
  }
  
  # Solving the optimization problem 
  solve(lp1)
  
  return(apply(matrix(get.variables(lp1)[1:(n * k)], ncol = k), 1, which.max))
  
}


# Examples ----------------------------------------------------------------


# Simulated data (n = 200)
{
  set.seed(2)
  
  # x: x-coordinate
  # y: y-coordinate
  # w: weights for objects (must be integers)
  # id: identification for each object
  koord <- data.frame(x = c(rnorm(50, mean = 1.5, sd = 0.8), rnorm(50, mean = 5.5, sd = 0.8), 
                            rnorm(50, mean = 3.0, sd = 0.8), rnorm(50, mean = 4.5, sd = 0.8)),
                      y = c(rnorm(50, mean = 1.0, sd = 0.8), rnorm(50, mean = 1.5, sd = 1.0), 
                            rnorm(50, mean = 4.5, sd = 1.0), rnorm(50, mean = 4.5, sd = 0.8)),
                      w =  round(runif(200, min = 0.51, max = 10.49)),
                      id = 1:200)
  n <- length(koord$x)
  plot(koord[, 1:2], cex = koord$w / 3, pch = 19, main = paste("Testidata, n =", n),
       xlim = c(min(koord$x) - 1, max(koord$x) + 1), ylim=c(min(koord$y) - 1 ,max(koord$y) + 1))
}

# Number of clusters
k <- 5

# Initial mu with k-means
{
  init_kmeans <- kmeans(cbind(rep(koord$x, koord$w), rep(koord$y, koord$w)), centers = k)
  
  init_mu <- init_kmeans$centers
  
  plot(koord[, 1:2], pch = 19, main = paste("Initial mu, n =", length(koord$x)), cex = koord$w / 3,
       xlim = c(min(koord$x) - 1, max(koord$x) + 1), ylim = c(min(koord$y) - 1, max(koord$y) + 1))
  points(init_mu[, 1], init_mu[, 2], cex = 5, pch = 4, lwd = 5, col = "red")
}

# Prior cluster sizes and corresponding probabilities
{
  # Mean
  pr_mean <- round(sum(koord$w) / k)
  
  # Standard deviation
  pr_sd <- 3
  
  # Max width for prior
  pr_width <- 30
  
  # All possible cluster sizes
  cl_size <- (pr_mean - pr_width):(pr_mean + pr_width)
  
  ## Uniform prior
  #d = rep(1/A,A)
  
  # Normal prior
  prob <- dnorm(cl_size, mean = pr_mean, sd = pr_sd)
}

# Function call
{
temp <- prob_clust(data = koord[, 1:2], weights = koord$w, k = k, init_mu = init_mu, 
                  prior_cl_sizes = cl_size, prior_prob = prob)
koord$cl <- temp[[1]]
mu <- temp[[2]]
}

# Picture 1
{
  plot(koord[, 1:2], cex = koord$w / 3, pch = 19, main = paste("Probabilistic Clustering, n =", n), col = c_col[koord$cl], 
       xlim = c(min(koord$x) - 2, max(koord$x) + 1), ylim = c(min(koord$y) - 1, max(koord$y) + 1))
  points(mu, cex = 2, pch = 4, lwd = 4)
  legend(min(koord$x) - 2, max(koord$y) + 1,col = c_col[1:k], pch = 19, cex = 1,
         legend = paste("Cluster", 1:k, paste("  (", apply(X = t(1:k), MARGIN = 2, FUN = function(x) {sum(koord$w[koord$cl == x])} ),")", sep="")))
}



