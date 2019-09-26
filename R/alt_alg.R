#' Full alternating algorithm
#'
#' @param coord The coordinates data.
#' @param weights Weights of the points in a vector.
#' @param k Number of clusters.
#' @param N Number of iterations.
#' @param limits Limits for the cluster size in a list.
#' @param d Distance function used in clustering.
#' @param lambda Outgroup parameter.
#' @param frac_memb Can points be partially allocated?
#' @param place_to_AP Place the cluster head in a point?
#' @param fixed_mu Possible fixed center locations.
#'
#' @return Clustering object with allocation, center locations and the value of the objective function
#' @export
#'
#' @examples
alt_alg <- function(coord, weights, k, N = 10, limits = bounds(data, k, radius = 100), d = euc_dist2, lambda = NULL, frac_memb = FALSE, place_to_AP = TRUE, fixed_mu = NULL){
  
  if(N < 2){
    N <- 2
  }
  
  w <- weights
  # Call prob_clust function
  temp <- prob_clust(data = coord,
                     weights = w,
                     k = k,
                     prior_dist = "uniform",
                     range = c(limits$L, limits$U),
                     d = d,
                     frac_memb = frac_memb,
                     fixed_mu = fixed_mu)
  # Print the number of completed laps
  print(paste("###################################################################"))
  print(paste("###################################################################"))
  print(paste("######################## Laps completed:", 1, "########################"))
  print(paste("###################################################################"))
  print(paste("###################################################################"))
  min_obj <- temp$obj
  best_temp <- temp
  for (i in 2:N) {
    temp <- prob_clust(data = coord,
                       weights = w,
                       k = k,
                       prior_dist = "uniform",
                       range = c(limits$L, limits$U),
                       d = d,
                       frac_memb = frac_memb,
                       fixed_mu = fixed_mu)
    
    # Print the number of completed laps
    print(paste("###################################################################"))
    print(paste("###################################################################"))
    print(paste("######################## Laps completed:", i, "########################"))
    print(paste("###################################################################"))
    print(paste("###################################################################"))
    if(temp$obj < min_obj){
      min_obj <- temp$obj
      best_temp <-  temp
    }
  }
  return(best_temp)
}