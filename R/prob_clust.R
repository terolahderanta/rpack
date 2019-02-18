#' @title Probabilistic Clustering
#'
#' @description Alternating algorithm for maximizing the joint density. Constraints on
#' cluster sizes are given as a prior distribution.
#'
#' @param data A matrix or data.frame containing the data.
#' @param weights A vector of weights for each data point.
#' @param init_mu Parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param prior_cl_sizes All the possible values for cluster sizes.
#' @param prior_prob The corresponding probabilities (sum to 1).
#' @param lambda Outgroup-parameter.
#' @return A list containting the new cluster allocations for each object in data, the new cluster center locations and maximum of the objective function.
#' @export prob_clust
prob_clust <- function(data, weights, k, init_mu, prior_cl_sizes, prior_prob, lambda = 0){

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

  max_sim <- 50
  for (iter in 1:max_sim) {
    # Old mu is saved to check for convergence
    old_mu <- mu

    # Clusters in equally weighted data (Allocation-step)
    temp_allocation <- prob_clust_estep(data_ew, mu, k, prior_cl_sizes, prior_prob, lambda)
    
    # Clusters in equally weighted data
    clusters_ew <- temp_allocation[[1]]
    
    # Outgroup-clusters (cluster 99)
    clusters_ew <- ifelse(clusters_ew == (k+1), 99, clusters_ew)
    
    # Maximum value of the objective function
    obj_max <- temp_allocation[[2]]

    # Updating cluster centers (Parameter-step)
    mu <- prob_clust_mstep(data_ew, clusters_ew, k)

    print(paste("Iteration:",iter))

    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }
  
  # Converting clusters_ew to original data
  for (i in 1:n) {
    clusters[i] <- getmode(clusters_ew[id_ew == i])
  }
  
  return(list(clusters, mu, obj_max))
}

#' Update the parameters (centers) for each cluster.
#'
#' @param data_ew A matrix or data.frame containing the data.
#' @param clusters_ew A vector of cluster assignments for each data point.
#' @param k The number of clusters.
#' @return New cluster centers.
prob_clust_mstep <- function(data_ew, clusters_ew, k){

  # Matrix for cluster centers
  mu <- matrix(0, nrow = k, ncol = length(data[1,]))

  # Update mu given
  for (i in 1:k) {
    mu[i,] <- c(mean(data_ew[clusters_ew == i, 1]), mean(data_ew[clusters_ew == i, 2]))
  }

  return(mu)
}

#' Update cluster allocations by maximizing the joint log-likelihood (Allocation-step).
#'
#' @param data_ew A matrix or data.frame containing the data, where each object is considered to be equally weighted.
#' @param mu Parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param prior_cl_sizes All the possible values for cluster sizes.
#' @param prior_prob Corresponding probabilities (sum to 1).
#' @param lambda Outgroup-parameter.
#' @return New cluster allocations for each object in data_ew
prob_clust_estep <- function(data_ew, mu, k, prior_cl_sizes, prior_prob, lambda){

  # Number of objects in data_ew
  n <- length(data_ew[,1])

  # Log priors for possible cluster sizes
  prior_log_prob <- log(prior_prob)

  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(mu, MARGIN = 1, FUN = mvtnorm::dmvnorm, x = data_ew, sigma = diag(2), log = TRUE)

  # Number of possible cluster sizes
  A <- length(prior_cl_sizes)

  # New linear program model object containing, where the number of decision variables is n * k + k * A
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = n * k + k * A)

  # Maximizing the joint log-likelihood
  lpSolveAPI::lp.control(lp1, sense = "max")

  # Binary integer problem
  lpSolveAPI::set.type(lp1, 1:(n*k + k*A), "binary")

  # Objective function in the optimization problem
  lpSolveAPI::set.objfn(lp1, c(C, t(matrix(1, ncol = k, nrow = A) * prior_log_prob)))

  # Constraint 1: each point is assigned to exactly one cluster
  for (i in 1:n) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(rep(1:n == i, k)), rep(0, k * A)), "=", 1)
  }

  # Constraint 2: each cluster has exactly one size
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, c(rep(0, k * n), as.numeric(rep(1:k == i, A))), "=", 1)
  }

  # Constraint 3: each cluster has S_a points
  temp <- c(t(matrix(1, ncol = n, nrow = k) * 1:k))
  temp2 <- c(t(-matrix(1, ncol = k, nrow = A) * prior_cl_sizes))
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(temp == i, n), temp2 * as.numeric(rep(1:k == i, A))), "=", 0)
  }

  # Solving the optimization problem
  solve(lp1)
  
  # Maximum of the objective function
  obj_max <- round(get.objective(lp1), digits = 2)

  # Print the value of the objective function
  print(paste("Value of the objective function:", obj_max))
  
  return(list(apply(matrix(lpSolveAPI::get.variables(lp1)[1:(n * k)], ncol = k), 1, which.max), obj_max))
}
