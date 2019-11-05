#' Full alternating algorithm
#'
#' @param coords Coordinates of the data points.
#' @param weights Weights of the points in a vector.
#' @param k Number of clusters.
#' @param N Number of iterations.
#' @param range Limits for the cluster size in a list.
#' @param capacity_weights Different weights for capacity limits.
#' @param d Distance function used in clustering.
#' @param mu_initialization Method to initialize mu, default is kmpp.
#' @param lambda Outgroup parameter.
#' @param frac_memb Can points be partially allocated?
#' @param place_to_point Place the cluster head in a point?
#' @param fixed_mu Possible fixed center locations.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param dist_mat Distance matrix for all the points. 
#' @param predet_locations Choose centers only from predetermined locations.
#'
#' @return Clustering object with allocation, center locations and the value of the objective function
#' @export
#'
#' @examples
alt_alg <- function(coords, weights, k, N = 10, range = as.numeric(bounds(weights, k, radius = 100)),
                    capacity_weights = weights, d = euc_dist2, mu_initialization = NULL, lambda = NULL,
                    frac_memb = FALSE, place_to_point = TRUE, fixed_mu = NULL, gurobi_params = NULL,
                    dist_mat = NULL, predet_locations = NULL){
  
  if(N < 2){
    N <- 2
  }
  
  w <- weights
  if(is.null(mu_initialization)){
    init_mu <- NULL
    
  } else if(mu_initialization == "random"){
    init_mu <- list()
    for (i in 1:N) {
      init_mu[[i]] <- as.matrix(coords[sample(x = 1:nrow(coords), size = k),])
    }
    
  } else {
    init_mu <- NULL
  }
  
  # Print the information about run
  cat(paste("Progress bar (N = ", N,"):\n", sep = ""))
  cat(paste("______________________________\n"))
  
  progress_bar <- 0
  
  # Call prob_clust function first time
  temp <- prob_clust(data = coords,
                     weights = w,
                     k = k,
                     prior_dist = "uniform",
                     range = c(range[1], range[2]),
                     capacity_weights = capacity_weights,
                     lambda = lambda,
                     d = d,
                     init_mu = init_mu[[1]],
                     place_to_point = place_to_point,
                     frac_memb = frac_memb,
                     fixed_mu = fixed_mu,
                     gurobi_params = gurobi_params,
                     dist_mat = dist_mat,
                     predet_locations = predet_locations)
  
  min_obj <- temp$obj
  best_temp <- temp
  for (i in 2:N) {
    temp <- prob_clust(data = coords,
                       weights = w,
                       k = k,
                       prior_dist = "uniform",
                       range = c(range[1], range[2]),
                       capacity_weights = capacity_weights,
                       lambda = lambda,
                       d = d,
                       init_mu = init_mu[[i]],
                       place_to_point = place_to_point,
                       frac_memb = frac_memb,
                       fixed_mu = fixed_mu,
                       gurobi_params = gurobi_params,
                       dist_mat = dist_mat,
                       predet_locations = predet_locations)
    
    # Print the number of completed laps
    if(floor((i/N)*30) > progress_bar){
      cat(paste0(rep("#", floor((i/N)*30) - progress_bar), collapse = ""))  
      progress_bar <- floor((i/N)*30)
    }
    
    # Save the iteration with the lowest value of objective function
    if(temp$obj < min_obj){
      min_obj <- temp$obj
      best_temp <-  temp
    }
  }
  cat("\n")
  return(best_temp)
}
