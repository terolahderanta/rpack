# Probabilistic clustering algorithm where constraints on cluster sizes are given as a prior distribution. 

# Packages ----------------------------------------------------------------
{
library(lpSolveAPI)
library(lpSolve)
library(mvtnorm)
library(pracma)
c_col = c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
          "chocolate","darkgrey","darkred","yellow3","darkgreen","wheat3","magenta",
          "palegreen2","violetred","seagreen2","tomato4","steelblue1","royalblue",
          "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
          "coral3","beige","brown4","indianred1","lightgreen","orchid")
}
# Functions ---------------------------------------------------------------
{
plot_clusters <- function(data, weights, clusters, mu, main){
  # Plots the clusters
  
  n <- length(data[,1])
  k <- length(mu[,1])
  
  plot(data, cex = weights / 3, pch = 19, main = paste(main, " (n = ", n, ")",sep = ""), col = c_col[clusters], 
       xlim = c(min(data[,1]) - 2, max(data[,1]) + 1), ylim = c(min(data[,2]) - 1, max(data[,2]) + 1))
  points(data[clusters == 99,])
  points(mu, cex = 2, pch = 4, lwd = 4)
  legend(min(data[,1]) - 2.5, max(data[,2]) +1,col = c_col[c(1:k, 99)], pch = 19, cex = 0.8,
         legend = paste("Cluster", c(1:k, 99), paste("  (", apply(X = t(c(1:k, 99)), MARGIN = 2, FUN = function(x) {sum(weights[clusters == x])} ),")", sep="")))
  
}

kmpp <- function(X, k) {
  #K-means++ algorithm
  
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  
  for (i in 2:k) {
    dm <- distmat(X, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  kmeans(X, X[C, ])
}

getmode <- function(v) {
  # Returns the mode of vector v
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

prob_clust <- function(data, weights, k, init_mu, prior_cl_sizes, prior_prob, lambda = 0){
  # Alternating algorithm for maximizing the joint density.
  # Parameters
  #   data: matrix or data.frame containing the data.
  #   weights: vector of weights for each data point.
  #   init_mu: parameters (locations) that define the k distributions.
  #   k: number of clusters.
  #   prior_cl_sizes: all the possible values for cluster sizes.
  #   prior_prob: corresponding probabilities (sum to 1).  
  #   lambda: outlier-parameter
  # Returns
  #   clusters: new cluster allocations for each object in data
  #   mu: new cluster center locations
  
  print("========== Step 1 ==========")
  
  # Number of objects in data
  n <- length(data[,1])
  
  # Weights must be integers
  weights <- round(weights)
  
  # Equally weighted data
  data_ew <- apply(data, MARGIN = 2, FUN = rep, times = weights)
  
  # Original id:s for each data point in equally weighted data
  id_ew <- rep(1:n, times = weights)
  
  # Initial clusters
  clusters <- rep(0,n)
  
  #Cluster centers
  mu <- init_mu
  
  # Maximum number of laps
  max_sim <- 30
  
  # Narrowing the prior interval
  narr <- round(max(weights)/2)
  
  # Narrowed prior probabilities
  prior_prob_narr <- c(rep(0,narr), prior_prob[(1+narr):(length(prior_prob)-narr)], rep(0,narr))
  
  for (iter in 1:max_sim) {
    # Old mu is saved to check for convergence
    old_mu <- mu
    
    # Clusters in equally weighted data (Allocation-step)
    temp_allocation <- prob_clust_allocation(data_ew, mu, k, prior_cl_sizes, prior_prob, lambda)
    clusters_ew <- temp_allocation[[1]]
    clusters_ew <- ifelse(clusters_ew == (k+1), 99, clusters_ew)
    obj_min <- temp_allocation[[2]]
    
    # Converting clusters_ew to original data
    #for (i in 1:n) {
    #  clusters[i] <- getmode(clusters_ew[id_ew == i])
    #}
    
    # Updating cluster centers (Parameter-step)
    mu <- prob_clust_parameter(data_ew, clusters_ew, k)
    
    print(paste("Iteration:",iter))
    
    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }
  
  # Converting clusters_ew to original data
  for (i in 1:n) {
    clusters[i] <- getmode(clusters_ew[id_ew == i])
  }
  
  print("========== Step 2 ==========")
  
  max_indiv <- 10
  
  for (indiv_lap in 1:max_indiv) {
    # Old mu is saved to check for convergence
    old_mu <- mu
    
    # Updating point allocation individually
    clusters <- prob_clust_allocation_indiv(data, weights, clusters, mu, k, prior_cl_sizes, prior_prob, lambda)
    
    # Updating cluster centers
    mu <- prob_clust_parameter_weights(data, clusters, weights, k)
    
    # Value of the objective function
    obj_min <- obj_function(data, weights, clusters, mu, prior_cl_sizes, prior_prob, lambda)
    
    print(paste("Indiv. update:", indiv_lap))
    print(paste("Objective function:", round(obj_min, digits = 2)))
    
    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }
  
  return(list(clusters, mu, obj_min))
}

prob_clust_allocation_indiv <- function(data, weights, clusters, mu, k, prior_cl_sizes, prior_prob, lambda){
  # Updates cluster allocations by individually allocationg points. 
  # Parameters
  #   data: matrix or data.frame containing the data, where each object is considered to be equally weighted.
  #   mu: parameters (locations) that define the k distributions.
  #   k: number of clusters.
  #   prior_cl_sizes: all the possible values for cluster sizes.
  #   prior_prob: corresponding probabilities (sum to 1).  
  # Returns
  #   new_allocation: new cluster allocations for each object in data_ew
  
  # Number of objects in data_ew
  n <- length(data[, 1])
  
  # Maximum number of update laps
  max_update <- 100
  
  # Cluster sizes initialization
  cluster_size = rep(0, k)
  
  # Lower limit for prior cluster size
  L <- prior_cl_sizes[1]
    
  # Upper limit for prior cluster size
  U <- prior_cl_sizes[length(prior_cl_sizes)]
  
  # Updating objects to clusters individually
  for(update_lap in 1:max_update){

    # To check convergence
    old_clusters <- clusters
    
    # Random order to allocate points
    points <- sample(1:n)
    
    # Matrix contains the log-likelihoods of the individual data points
    C <- apply(mu, MARGIN = 1, FUN = dmvnorm, x = data, sigma = diag(2), log = TRUE)
    
    # One allocation update for all the points 
    for (i in points) {
      
      # Cluster sizes in the current allocation
      for (j in 1:k) {
        cluster_size[j] <- sum(weights[clusters == j])  
      }
      
      # If the current point is an outlier, ignore it
      if(clusters[i] != 99) {
        
        # Check to see, if the point can be removed from the original cluster
        if((cluster_size[clusters[i]] - weights[i]) >= L) {
          
          # Original cluster
          orig_cluster = clusters[i]
          
          temp_clust <- matrix(clusters, nrow = n, ncol = k)
          temp_clust[i,] <- 1:k
          
          # Check to see, which clusters have enough space for the current point
          acc_move <- ((cluster_size + weights[i]) <= U)
          
          # The point can be allocated to the original cluster
          acc_move[orig_cluster] <- TRUE
          
          # Calculate the value of the objective function 
          temp_clust <- temp_clust[,acc_move]
          temp_obj <- apply(X = temp_clust, MARGIN = 2, FUN = obj_function, 
                            data = data, weights = weights, mu = mu,
                            prior_cl_sizes = prior_cl_sizes, 
                            prior_prob = prior_prob, lambda = lambda)
          
          # Choose the highest value
          clusters[i] <- (1:k)[acc_move][which.max(temp_obj)]
        }
      }
    }
    # If nothing is changing, stop
    if(all(old_clusters == clusters)) break
    
  }
  return(clusters)
  
}

obj_function <- function(data, weights, clusters, mu, prior_cl_sizes, prior_prob, lambda){
  # Returns the objective function value given clusters and mu
  # Parameters
  #   data: matrix or data.frame containing the data.
  #   clusters: current point allocation to clusters
  #   mu: parameters (locations) that define the k distributions.
  #   prior_cl_sizes: all the possible values for cluster sizes.
  #   prior_prob: corresponding probabilities (sum to 1).  
  # Returns
  #   obj: value of the objective function.
  
  n <- length(data[,1])
  k <- length(mu[,1])
  
  # Log priors for possible cluster sizes
  prior_log_prob <- log(prior_prob)
  
  # Lower limit for prior cluster size
  L <- prior_cl_sizes[1]
  
  # Upper limit for prior cluster size
  U <- prior_cl_sizes[length(prior_cl_sizes)]
  
  # Cluster sizes initialization
  cluster_size = rep(0, k)
  
  # Cluster sizes in the current allocation
  for (j in 1:k) {
    cluster_size[j] <- sum(weights[clusters == j])  
  }
  
  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(mu, MARGIN = 1, FUN = dmvnorm, x = data, sigma = diag(2), log = TRUE)
  
  # Scaling the tuning parameter lambda
  nu2 <- -(mean(dist(data))/sqrt(k)) ^ 2
  
  # Cluster allocation
  z <- matrix(0, nrow = n, ncol = k)
  
  for (i in 1:k) {
    z[,i] <- as.numeric(clusters == i) 
  }
  
  # Cluster size allocation
  B <- matrix(0, nrow = k, ncol = length(prior_cl_sizes))
  
  for (i in 1:length(prior_cl_sizes)) {
    B[,i] <- as.numeric(cluster_size == prior_cl_sizes[i]) 
  }
  
  # Outlier allocation
  z_out <- as.numeric(clusters == 99)
  
  Cz = c(C * z)
  
  # Objective function
  return(sum(Cz * weights) + sum(prior_log_prob * t(B)) + nu2*lambda*sum(z_out * weights))
}

prob_clust_parameter_weights <- function(data, clusters, weights, k){
  # Updates the parameters (centers) for each cluster.
  # Parameters
  #   data: matrix or data.frame containing the data.
  #   clusters: vector of cluster assignments for each data point.
  #   weights: weights of the data points 
  #   k: number of clusters.
  # Returns
  #   mu: new cluster centers.
  
  # Matrix for cluster centers
  mu <- matrix(0,nrow = k, ncol = length(data[1,]))
  
  # Update mu for each cluster
  for (i in 1:k) {
    # Weighted mean
    mu[i,] <- apply(data[clusters == i,]*weights[clusters == i], 2, FUN=sum) / sum(weights[clusters == i])  
  }
  
  return(mu)
}

prob_clust_parameter <- function(data_ew, clusters_ew, k){
  # Updates the parameters (centers) for each cluster.
  # Parameters
  #   data_ew: matrix or data.frame containing the data.
  #   clusters_ew: vector of cluster assignments for each data point.
  #   k: number of clusters.
  # Returns
  #   mu: new cluster centers.
  
  # Matrix for cluster centers
  mu <- matrix(0,nrow = k, ncol = length(data_ew[1,]))
  
  # Update mu
  for (i in 1:k) {
    mu[i,] <- c(mean(data_ew[clusters_ew == i,1]),mean(data_ew[clusters_ew == i,2]))
    #mu[i,] <- apply(data[clusters == i,] * weights[clusters == i], 2, FUN=sum) / sum(weights[clusters == i])  
  }
  
  return(mu)
}

#TODO: make lambda do the outlier thingy
prob_clust_allocation <- function(data_ew, mu, k, prior_cl_sizes, prior_prob, lambda){
  # Updates cluster allocations by maximizing the joint log-likelihood (M-step).
  # Parameters
  #   data_ew: matrix or data.frame containing the data, where each object is considered to be equally weighted.
  #   mu: parameters (locations) that define the k distributions.
  #   k: number of clusters.
  #   prior_cl_sizes: all the possible values for cluster sizes.
  #   prior_prob: corresponding probabilities (sum to 1).
  #   lambda: outlier-parameter
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
  
  obj_min <- round(get.objective(lp1), digits = 2)
  
  # Print the value of the objective function
  print(paste("Value of the objective function:", obj_min))
  
  return(list(apply(matrix(get.variables(lp1)[1:(n * k)], ncol = k), 1, which.max), obj_min))
  
}
}

# Examples ----------------------------------------------------------------


# Simulated data (n = 200)
{
  set.seed(4)
  
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
k <- 10

# Initial mu with k-means
{
  init_kmpp <- kmpp(cbind(rep(koord$x, koord$w), rep(koord$y, koord$w)), k)
  init_mu <- init_kmpp$centers
  
  plot(koord[, 1:2], pch = 19, main = paste("Initial mu, n =", length(koord$x)), cex = koord$w / 3,
       xlim = c(min(koord$x) - 1, max(koord$x) + 1), ylim = c(min(koord$y) - 1, max(koord$y) + 1))
  points(init_mu[, 1], init_mu[, 2], cex = 5, pch = 4, lwd = 5, col = "red")
}

# Prior cluster sizes and corresponding probabilities
{
  # Mean
  pr_mean <- round(sum(koord$w) / k)
  
  # Standard deviation
  pr_sd <- 5
  
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
  plot_clusters(koord[,1:2], koord$w, koord$cl, mu, main = paste("Range: ", min(cl_size), "-", max(cl_size), sep = ""))
}

