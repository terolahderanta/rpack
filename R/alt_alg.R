#' Full alternating algorithm
#'
#' @param coords Coordinates of the data points.
#' @param weights Weights of the points in a vector.
#' @param k Number of clusters.
#' @param N Number of iterations.
#' @param range Limits for the cluster size in a list.
#' @param capacity_weights Different weights for capacity limits.
#' @param d Distance function used in clustering.
#' @param center_init Options to initialize center locations. Default is "random" and other choice is "kmpp". 
#' @param lambda Outgroup parameter.
#' @param frac_memb Can points be partially allocated?
#' @param place_to_point Place the cluster head in a point?
#' @param fixed_centers Possible fixed center locations.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param multip_centers Vector (n-length) defining how many centers a point is allocated to.
#' @param parallel Logical indicator to use parallel computing.
#' @param dist_mat Distance matrix for all the points. 
#' @param print_output Different types of printing outputs, 1 is default, 2 stepwise-print and 3 is full print.
#'
#' @return Clustering object with allocation, center locations and the value of the objective function
#' @export
#'
#' @examples
alt_alg <- function(coords, 
                     weights, 
                     k, 
                     N = 10, 
                     range,
                     capacity_weights = weights, 
                     d = euc_dist2, 
                     center_init = "random", 
                     lambda = NULL,
                     frac_memb = FALSE, 
                     place_to_point = TRUE, 
                     fixed_centers = NULL, 
                     gurobi_params = NULL,
                     multip_centers = rep(1, nrow(coords)),
                     parallel = FALSE,
                     dist_mat = NULL,
                     print_output = 1){
  
  # Print the information about run
  if(print_output == 1){
    cat(paste("Progress (N = ", N,"):\n", sep = ""))
    cat(paste("______________________________\n"))
    progress_bar <- 0
  }
  
  # Calculate distance matrix
  if(is.null(dist_mat) & place_to_point){
    
    # Calculate distances with distance metric d
    dist_mat <- apply(
      X = coords,
      MARGIN = 1,
      FUN = function(x) {
        apply(
          X = coords,
          MARGIN = 1,
          FUN = d,
          x2 = x
        )
      }
    )
    
    # Normalizing distances
    dist_mat <- dist_mat / max(dist_mat)
  } else {
    dist_mat <- NULL
  }
  
  # Normalization for the capacity weights
  max_cap_w <- max(capacity_weights)
  capacity_weights <- capacity_weights/max_cap_w
  range <- range/max_cap_w
  
  # Normalization for the weights
  weights <- weights/max(weights)
  
  
  for (i in 1:N) {
    
    if(print_output == 2){
      cat(paste("\nIteration ", i, "\n----------\n", sep = ""))
    }
    
    # One clustering
    temp_clust <- capacitated_LA(coords = coords,
                                 weights = weights,
                                 k = k,
                                 ranges = range,
                                 capacity_weights = capacity_weights,
                                 lambda = lambda,
                                 d = d,
                                 dist_mat = dist_mat,
                                 center_init = center_init,
                                 place_to_point = place_to_point,
                                 frac_memb = frac_memb,
                                 fixed_centers = fixed_centers,
                                 gurobi_params = gurobi_params,
                                 multip_centers = multip_centers,
                                 parallel = parallel)
    
    # Save the first iteration as the best one
    if(i == 1){
      min_obj <- temp_clust$obj
      best_clust <- temp_clust
    }
    
    # Print the number of completed laps
    if(print_output == 1 & (floor((i/N)*30) > progress_bar)){
      cat(paste0(rep("#", floor((i/N)*30) - progress_bar), collapse = ""))  
      progress_bar <- floor((i/N)*30)
    } 
    
    # Save the iteration with the lowest value of objective function
    if(temp_clust$obj < min_obj){
      min_obj <- temp_clust$obj
      best_clust <-  temp_clust
    }
  }
  
  cat("\n")
  return(best_clust)
}