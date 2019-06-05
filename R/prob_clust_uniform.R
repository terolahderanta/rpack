#' Probabilistic Clustering with Simple Parameters
#'
#' Alternating algorithm for maximizing the joint density. This is a simpler version of
#' the probabilistic clustering algorithm prob_clust, with the constraints to cluster sizes given
#' only as a boundary from L to U.
#'
#' @param data A matrix or data.frame containing the data.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param init_mu Parameters (locations) that define the k distributions.
#' @param L The lower limit for cluster sizes.
#' @param U The upper limit for cluster sizes.
#' @param fixed_mu Predetermined center locations.
#' @param lambda Outgroup-parameter.
#' @param d The distance function.
#' @return A list containing cluster allocation, cluster center and the current value of the objective function.
#' @keywords internal
prob_clust_uniform <- function(data, weights, k, init_mu, L, U, fixed_mu = NULL, lambda = NULL, d = euc_dist2){

  # Number of objects in data
  n <- nrow(data)
  
  # Weights must be integers
  weights <- round(weights)
  
  # Initial clusters
  clusters <- rep(0,n)
  
  # Number of fixed centers
  n_fixed <- ifelse(is.null(fixed_mu), 0, nrow(fixed_mu))
  
  if(n_fixed > 0){
    # Cluster centers with fixed centers first
    mu <- rbind(fixed_mu, init_mu[(n_fixed + 1):k,])
  } else {
    # Cluster centers without any fixed centers
    mu <- init_mu
  }
  # Maximum number of laps
  max_sim <- 50
  
  for (iter in 1:max_sim) {
    # Old mu is saved to check for convergence
    old_mu <- mu
    
    # Clusters in equally weighted data (Allocation-step)
    temp_allocation <- prob_clust_allocation_weights(data, weights, mu, k, L, U, lambda, d = d)
    assign_frac <- temp_allocation[[1]]
    obj_max <- temp_allocation[[2]]
    
    # Updating cluster centers (Parameter-step)
    mu <- prob_clust_parameter_weights(data, assign_frac, weights, k, fixed_mu = fixed_mu, d = d)
    
    print(paste("Iteration:",iter))
    
    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }
  
  # Hard clusters from assign_frac
  clusters <- apply(assign_frac, 1, which.max)
  
  # Cluster 99 is the outgroup
  clusters <- ifelse(clusters == (k+1), 99, clusters)
  
  # Return cluster allocation, cluster center and the current value of the objective function
  return(list(clusters = clusters, centers = mu, obj = obj_max, assign_frac = assign_frac))
}

#' Calculate the objective function value, given clusters and mu.
#'
#' @param data A matrix or data.frame containing the data.
#' @param weights A vector of weights for each data point.
#' @param clusters Current point allocation to clusters
#' @param mu Parameters (locations) that define the k distributions.
#' @param lambda Outgroup-parameter.
#' @return The value of the objective function.
#' @keywords internal
obj_function <- function(data, weights, clusters, mu, lambda){

  n <- nrow(data)
  k <- nrow(mu)

  # Matrix contains the log-likelihoods of the individual data points
  C <- apply(
    mu,
    MARGIN = 1,
    FUN = mvtnorm::dmvnorm,
    x = data,
    sigma = diag(2),
    log = TRUE
  )

  # Scaling the tuning parameter lambda
  nu2 <- -(mean(stats::dist(data))/sqrt(k)) ^ 2

  # Cluster allocation
  z <-  sapply(
    1:k,
    FUN = function(x) as.numeric(clusters == x)
  )

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
#' @keywords internal
prob_clust_parameter <- function(data_ew, clusters_ew, k){

  # Matrix for cluster centers
  mu <- matrix(0, nrow = k, ncol = ncol(data_ew))

  #col_means <- function(data) {
  #  colSums(data) / nrow(data)
  #}

  #sapply(1:k, function(i) {
  #  mu[i, ] = col_means(data_ew[clusters_ew == i])
  #})

  #mu <- apply(data_ew, 2, function(data_col) {
  #  mean(data_col[clusters_ew == i])
  #})

  # Update mu
  for (i in 1:k) {
    mu[i,] <- c(mean(data_ew[clusters_ew == i, 1]),
                mean(data_ew[clusters_ew == i, 2]))
  }

  return(mu)
}

#' Update cluster allocations by maximizing the joint log-likelihood.
#'
#' @param data A matrix or data.frame containing the data, where each object is considered to be equally weighted.
#' @param weights A vector of weights for each data point.
#' @param mu The parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param L The lower limit for cluster sizes.
#' @param U The upper limit for cluster sizes.
#' @param lambda Outgroup-parameter
#' @param d The distance function.
#' @return New cluster allocations for each object in data and the maximum of the objective function.
#' @keywords internal
prob_clust_allocation_weights <- function(data, weights, mu, k, L, U, lambda, d = euc_dist2){
  
  # Number of objects in data
  n <- nrow(data)
  
  # Is there an outgroup cluster
  is_outgroup <- !is.null(lambda)
  
  # Number of decision variables
  n_decision <- ifelse(is_outgroup, n * k + n, n * k)
  
  # Matrix contains the log-likelihoods of the individual data points
  C <- matrix(0, ncol = k, nrow = n)
  # TODO: Ways to measure distance with matrices
  for(i in 1:k){
    C[,i] <- apply(data, MARGIN = 1, FUN = d, x2 = mu[i,])
  }
  
  # New linear program model object
  lp1 <- lpSolveAPI::make.lp(nrow = 0, ncol = n_decision)
  
  # Maximizing the joint log-likelihood
  lpSolveAPI::lp.control(lp1, sense = "min") # FIXME: min or max depending on the distance measure 
  
  # Binary integer problem
  # TODO: This defines if we use hard or fractional clustering
  lpSolveAPI::set.type(lp1, 1:n_decision, "real")
  
  # Scaling the tuning parameter lambda
  nu2 <- -(mean(stats::dist(data))/sqrt(k)) ^ 2
  
  # Objective function in the optimization problem
  lpSolveAPI::set.objfn(lp1, c(C * weights, rep(nu2 * lambda, n) * weights))
  
  # Constraint 1: each point is assigned to exactly one cluster
  for (i in 1:n) {
    lpSolveAPI::add.constraint(lp1, rep(1,ifelse(is_outgroup, k + 1, k)), "=", 1,
                               indices = seq(from = i, to = n*(k-1) + i, by = n))
  }
  
  # Constraint 2: each cluster size must be higher than L
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, rep(1, n)*weights, ">=", L,
                               indices = ((i - 1) * n + 1):(n * i))
  }
  
  # Constraint 3: each cluster size must be lower than U
  for (i in 1:k) {
    lpSolveAPI::add.constraint(lp1, rep(1, n)*weights, "<=", U,
                               indices = ((i - 1) * n + 1):(n * i))
  }
  
  # Constraint 4: Each y_nk must be greater than or equal to 0
  for (i in 1:k) {
    for (j in 1:n) {
      lpSolveAPI::add.constraint(lp1, 1, ">=", 0,
                                 indices = (i - 1) * n + j)
    }
  }
  
  # Solving the optimization problem
  lpSolveAPI::solve(lp1)
  
  # Maximum/minimum of the objective function
  obj_max <- round(lpSolveAPI::get.objective(lp1), digits = 3)
  
  # Print the value of the objective function
  print(paste("Value of the objective function:", obj_max))
  
  return(list(matrix(lpSolveAPI::get.variables(lp1)[1:n_decision], ncol = ifelse(is_outgroup, k + 1, k)),
              obj_max))
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
#' @keywords internal
prob_clust_allocation <- function(data_ew, mu, k, L, U, lambda){

  # Number of objects in data_ew
  n <- nrow(data_ew)

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
#' @param weights The weigths of the objects in data.
#' @param clusters A vector of cluster assignments for each data point.
#' @param mu The parameters (locations) that define the k distributions.
#' @param k The number of clusters.
#' @param L A lower limit for cluster sizes.
#' @param U An upper limit for cluster sizes.
#' @param lambda Outgroup-parameter
#' @return New cluster allocations for each object in data_ew.
#' @keywords internal
prob_clust_allocation_indiv <- function(data, weights, clusters, mu, k, L, U, lambda){

  # Number of objects in data_ew
  #n <- length(data[,1])
  n <- nrow(data)

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
    C <- apply(
      mu,
      MARGIN = 1,
      FUN = mvtnorm::dmvnorm,
      x = data,
      sigma = diag(2),
      log = TRUE
    )

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
        temp_obj <- apply(
          X = temp_clust,
          MARGIN = 2,
          FUN = obj_function,
          data = data,
          weights = weights,
          mu = mu,
          lambda = lambda
        )

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
#' @param data A data.frame containing the data points, one per row.
#' @param clusters A vector of cluster assignments for each data point.
#' @param weights The weights of the data points
#' @param k The number of clusters.
#' @param fixed_mu Predetermined center locations.
#' @param d The distance function.
#' @return New cluster centers.
#' @keywords internal
prob_clust_parameter_weights <- function(data, clusters, weights, k, fixed_mu = NULL, d = euc_dist2){
  
  # Matrix for cluster centers
  mu <- matrix(0, nrow = k, ncol = ncol(data))
  
  # Number of fixed centers
  n_fixed <- ifelse(is.null(fixed_mu), 0, nrow(fixed_mu))
  
  if(n_fixed > 0){
    # Insert the fixed mu first
    for(i in 1:n_fixed){
      mu[i,] <- fixed_mu[i,]
    }
  }
  
  # Update mu for each cluster
  for (i in (ifelse(n_fixed > 0, n_fixed + 1, 1)):k) {
    # Weighted mean
    #mu[i,] <- colSums(data * weights * clusters[, i]) / sum(clusters[, i]*weights)
    
    # Compute medoids only with points that are relevant in the cluster i
    relevant_cl <- clusters[, i] > 0.001
    
    # Computing medoids for cluster i
    #   New weights are combined from assignment fractionals and weights
    mu[i,] <- as.matrix(medoid(data = data[relevant_cl,],
                               w = clusters[relevant_cl, i] * weights[relevant_cl],
                               d = d))
  }
  
  
  return(mu)
}
