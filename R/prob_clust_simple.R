#' @title Probabilistic Clustering with Simple Parameters
#'
#' @description Alternating algorithm for maximizing the joint density. This is a simpler version of
#' the probabilistic clustering algorithm prob_clust, with the constraints to cluster sizes given
#' only as a boundary from L to U.
#'
#' @param data A matrix or data.frame containing the data.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param init_mu Parameters (locations) that define the k distributions.
#' @param L The lower limit for cluster sizes.
#' @param U The upper limit for cluster sizes.
#' @param lambda Outgroup-parameter.
#' @return A list containing cluster allocation, cluster center and the current value of the objective function.
#' @export prob_clust_simple
prob_clust_simple <- function(data, weights, k, init_mu, L, U, lambda = NULL){

  print("========== Step 1 ==========")

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

  # Maximum number of laps
  max_sim <- 30

  # Narrowing the prior interval
  narr <- round(max(weights)/2)

  for (iter in 1:max_sim) {
    # Old mu is saved to check for convergence
    old_mu <- mu

    # Clusters in equally weighted data (Allocation-step)
    temp_allocation <- prob_clust_allocation(data_ew, mu, k, L + narr, U - narr, lambda)
    clusters_ew <- temp_allocation[[1]]
    clusters_ew <- ifelse(clusters_ew == (k+1), 99, clusters_ew)
    obj_max <- temp_allocation[[2]]

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
    clusters <- prob_clust_allocation_indiv(data, weights, clusters, mu, k, L, U, lambda)

    # Updating cluster centers
    mu <- prob_clust_parameter_weights(data, clusters, weights, k)

    # Value of the objective function
    obj_max <- obj_function(data, weights, clusters, mu, lambda)

    print(paste("Indiv. update:", indiv_lap))
    print(paste("Objective function:", round(obj_max, digits = 2)))

    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }

  # Return cluster allocation, cluster center and the current value of the objective function
  return(list(clusters, mu, obj_max))
}

#' Calculate the objective function value, given clusters and mu.
#'
#' @param data A matrix or data.frame containing the data.
#' @param weights FIXME
#' @param clusters Current point allocation to clusters
#' @param mu Parameters (locations) that define the k distributions.
#' @param lambda FIXME
#' @return The value of the objective function.
obj_function <- function(data, weights, clusters, mu, lambda){

  n <- length(data[,1])
  k <- length(mu[,1])

  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(mu, MARGIN = 1, FUN = mvtnorm::dmvnorm, x = data, sigma = diag(2), log = TRUE)

  # Scaling the tuning parameter lambda
  nu2 <- -(mean(stats::dist(data))/sqrt(k)) ^ 2

  # Cluster allocation
  z <- matrix(0, nrow = n, ncol = k)

  for (i in 1:k) {
    z[,i] <- as.numeric(clusters == i)
  }

  # Outlier allocation
  z_out <- as.numeric(clusters == 99)

  Cz = c(C * z)

  # Objective function
  return(sum(Cz * weights) + ifelse(!is.null(lambda), nu2 * lambda  *sum(z_out * weights), 0))
}

#' Update the parameters (centers) for each cluster.
#'
#' @param data_ew A matrix or data.frame containing the data.
#' @param clusters_ew A vector of cluster assignments for each data point.
#' @param k The number of clusters.
#' @return New cluster centers.
prob_clust_parameter <- function(data_ew, clusters_ew, k){

  # Matrix for cluster centers
  mu <- matrix(0,nrow = k, ncol = length(data_ew[1,]))

  # Update mu
  for (i in 1:k) {
    mu[i,] <- c(mean(data_ew[clusters_ew == i,1]),mean(data_ew[clusters_ew == i,2]))
    #mu[i,] <- apply(data_ew[clusters_ew == i,], 2, FUN=sum) / length(data_ew[clusters_ew == i,])
  }

  return(mu)
}

#' Update cluster allocations by maximizing the joint log-likelihood.
#'
#' @param data_ew A matrix or data.frame containing the data, where each object is considered to be equally weighted.
#' @param mu The parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param L The lower limit for cluster sizes.
#' @param U The upper limit for cluster sizes.
#' @param lambda Outgroup-parameter
#' @return New cluster allocations for each object in data_ew and the maximum of the objective function.
prob_clust_allocation <- function(data_ew, mu, k, L, U, lambda){

  # Number of objects in data_ew
  n <- length(data_ew[,1])
  
  # Is there an outgroup cluster
  is_outgroup <- !is.null(lambda)
  
  # Number of decision variables
  n_decision <- ifelse(is_outgroup, n * k + n, n * k)

  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(mu, MARGIN = 1, FUN = mvtnorm::dmvnorm, x = data_ew, sigma = diag(2), log = TRUE)

  # New linear program model object
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = n_decision)

  # Maximizing the joint log-likelihood
  lpSolveAPI::lp.control(lp1, sense = "max")

  # Binary integer problem
  lpSolveAPI::set.type(lp1, 1:n_decision, "binary")

  # Scaling the tuning parameter lambda
  nu2 <- -(mean(stats::dist(data_ew))/sqrt(k)) ^ 2

  # Objective function in the optimization problem
  lpSolveAPI::set.objfn(lp1, c(C, rep(nu2 * lambda, n)))

  # Constraint 1: each point is assigned to exactly one cluster
  for (i in 1:n) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(rep(1:n == i, ifelse(is_outgroup, k + 1, k)))), "=", 1)
  }

  # Constraint 2: each cluster size must be higher than L
  temp <- c(t(matrix(1, ncol = n, nrow = k) * 1:k))
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(temp == i, n), rep(0, ifelse(is_outgroup, n, 0))), ">=", L)
  }

  # Constraint 3: each cluster size must be lower than U
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(temp == i, n), rep(0, ifelse(is_outgroup, n, 0))), "<=", U)
  }

  # Solving the optimization problem
  solve(lp1)

  obj_max <- round(lpSolveAPI::get.objective(lp1), digits = 2)

  # Print the value of the objective function
  print(paste("Value of the objective function:", obj_max))

  return(list(apply(matrix(lpSolveAPI::get.variables(lp1)[1:n_decision], ncol = ifelse(is_outgroup, k + 1, k)), 1, which.max), obj_max))

}

#' Updates cluster allocations by individually allocationg points.
#'
#' @param data A matrix or data.frame containing the data, where each object is considered to be equally weighted.
#' @param weights FIXME
#' @param clusters FIXME
#' @param mu The parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param L A lower limit for cluster sizes.
#' @param U An upper limit for cluster sizes.
#' @param lambda Outgroup-parameter
#' @return New cluster allocations for each object in data_ew.
prob_clust_allocation_indiv <- function(data, weights, clusters, mu, k, L, U, lambda){

  # Number of objects in data_ew
  n <- length(data[,1])

  # Maximum number of update laps
  max_update <- 100

  # Cluster sizes initialization
  cluster_size = rep(0, k)

  # Updating objects to clusters individually
  for(update_lap in 1:max_update){
    #print(paste("Iteration:",update_lap))
    #print(paste("Objective function:",obj_function(data, weights, clusters, mu, lambda)))
    # To check convergence
    old_clusters <- clusters

    # Random order to allocate points
    points <- sample(1:n)

    # Matrix contains the log-likelihoods of the individual data points
    C <- apply(mu, MARGIN = 1, FUN = mvtnorm::dmvnorm, x = data, sigma = diag(2), log = TRUE)

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
                             data = data, weights = weights, mu = mu, lambda = lambda)

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

#' Updates the parameters (centers) for each cluster.
#'
#' @param data A matrix or data.frame containing the data.
#' @param clusters A vector of cluster assignments for each data point.
#' @param weights The weights of the data points
#' @param k The number of clusters.
#' @return New cluster centers.
prob_clust_parameter_weights <- function(data, clusters, weights, k){

  # Matrix for cluster centers
  mu <- matrix(0,nrow = k, ncol = length(data[1,]))

  # Update mu for each cluster
  for (i in 1:k) {
    # Weighted mean
    mu[i,] <- apply(data[clusters == i,]*weights[clusters == i], 2, FUN=sum) / sum(weights[clusters == i])
  }

  return(mu)
}

#' Updates cluster allocations by maximizing the joint log-likelihood (M-step).
#'
#' @param data A matrix or data.frame containing the data, where each object is considered to be equally weighted.
#' @param weights The weigths of the objects in data.
#' @param mu The parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param L A lower limit for cluster sizes.
#' @param U An upper limit for cluster sizes.
#' @param lambda FIXME
#' @return New cluster allocations for each object in data_ew
prob_clust_allocation_weights <- function(data, weights, mu, k, L, U, lambda){

  # Number of objects in data_ew
  n <- length(data[,1])

  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(mu, MARGIN = 1, FUN = mvtnorm::dmvnorm, x = data, sigma = diag(2), log = TRUE)

  # New linear program model object containing, where the number of decision variables is n * k + n
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = n * k + n)

  # Maximizing the joint log-likelihood
  lpSolveAPI::lp.control(lp1, sense = "max")

  # Binary integer problem
  lpSolveAPI::set.type(lp1, 1:(n * k + n), "binary")

  # Scaling the tuning parameter lambda
  nu2 <- -(mean(stats::dist(data))/sqrt(k)) ^ 2

  # Objective function in the optimization problem
  lpSolveAPI::set.objfn(lp1, c(C * weights, rep(nu2 * lambda, n)))

  # Constraint 1: each point is assigned to exactly one cluster
  for (i in 1:n) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(rep(1:n == i, k + 1))), "=", 1)
  }

  # Constraint 2: each cluster size must be higher than L
  temp <- c(t(matrix(1, ncol = n, nrow = k) * 1:k))
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(temp == i, n) * weights, rep(0, n)), ">=", L)
  }

  # Constraint 3: each cluster size must be lower than U
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, c(as.numeric(temp == i, n) * weights, rep(0, n)), "<=", U)
  }

  # Solving the optimization problem
  solve(lp1)

  obj_min <- round(lpSolveAPI::get.objective(lp1), digits = 2)

  # Print the value of the objective function
  print(paste("Value of the objective function:", obj_min))

  return(list(apply(matrix(lpSolveAPI::get.variables(lp1)[1:(n * k + n)], ncol = k + 1), 1, which.max), obj_min))

}



