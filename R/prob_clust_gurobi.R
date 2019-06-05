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
#' @param d The distance function.
#' @param fixed_mu Predetermined center locations.
#' @param lambda Outgroup-parameter.
#' @return A list containing cluster allocation, cluster center and the current value of the objective function.
#' @keywords internal
prob_clust_gurobi <- function(data, weights, k, init_mu, L, U, d = euc_dist2, fixed_mu = NULL, lambda = NULL){
  
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
    temp_allocation <- prob_clust_allocation_weights_gurobi(data, weights, mu, k, L, U, lambda, d = d)
    assign_frac <- temp_allocation[[1]]
    obj_max <- temp_allocation[[2]]
    
    # Updating cluster centers (Parameter-step)
    mu <- prob_clust_parameter_weights_gurobi(data, assign_frac, weights, k, fixed_mu, d = d)
    
    print(paste("Iteration:",iter))
    
    # If nothing is changing, stop
    if(all(old_mu == mu)) break
  }
  
  # Hard clusters from assign_frac
  clusters <- apply(assign_frac, 1, which.max)
  
  # Cluster 99 is the outgroup
  clusters <- ifelse(clusters == (k+1), 99, clusters)
  
  # Return cluster allocation, cluster center and the current value of the objective function
  return(list(clusters = clusters, centers = mu, obj = obj_max))
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
#' @return New cluster allocations for each object in data and the maximum of the objective function.
#' @keywords internal
prob_clust_allocation_weights_gurobi <- function(data, weights, mu, k, L, U, lambda, d = euc_dist2){
  
  # Number of objects in data
  n <- nrow(data)
  
  # Is there an outgroup cluster
  is_outgroup <- !is.null(lambda)
  
  # Number of decision variables
  n_decision <- ifelse(is_outgroup, n * k + n, n * k)
  
  # Matrix contains the log-likelihoods of the individual data points
  #C <- apply(mu, MARGIN = 1, FUN = mvtnorm::dmvnorm, x = data, sigma = diag(2), log = TRUE)
  
  C <- matrix(0, ncol = k, nrow = n)
  # TODO: Tähän tilalle muita keinoja mitata etäisyyttä
  for(i in 1:k){
    C[,i] <- apply(data, MARGIN = 1, FUN = d, x2 = mu[i,])
  }
  
  # First constraint
  const1 <- NULL
  for (i in 1:n) {
    const1 <- c(const1, as.numeric(rep(1:n == i, ifelse(is_outgroup, k + 1, k))))
  }
  
  # Second constraint
  const2 <- NULL
  temp <- c(t(matrix(1, ncol = n, nrow = k) * 1:k))
  for (i in 1:k) {
    const2 <- c(const2, as.numeric(temp == i, n)*weights, rep(0, ifelse(is_outgroup, n, 0)))
  }
  
  # Third constraint
  const3 <- const2
  
  # Gurobi-model
  model <- list()
  
  # Constraint matrix A
  # TODO: Use sparse matrix (function spMatrix)!!!!
  model$A          <- matrix(c(const1, const2, const3), 
                             nrow=n + 2*k, 
                             ncol=n_decision, 
                             byrow=T)
  
  model$obj        <- c(C * weights)
  
  model$modelsense <- 'min'
  
  model$rhs        <- c(rep(1, n), 
                        rep(L, k), 
                        rep(U, k))
  
  model$sense      <- c(rep('=', n), 
                        rep('>', k), 
                        rep('<', k))
  
  model$vtype      <- 'B'
  
  #params <- list(OutputFlag=0)
  
  # Solving the linear program
  result <- gurobi(model, params = NULL)
  
  assign_frac <- matrix(print(result$x), ncol = ifelse(is_outgroup, k + 1, k))
  
  obj_max <- round(result$objval, digits = 3)
  
  
  # Print the value of the objective function
  print(paste("Value of the objective function:", obj_max))
  
  # Clear space
  rm(model, result)
  
  return(list(assign_frac,
              obj_max))
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
prob_clust_parameter_weights_gurobi <- function(data, clusters, weights, k, fixed_mu = NULL, d = euc_dist2){
  
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
