#' Function to simulate coordinate and weights data
#'
#' @param n Number of points.
#' @param w_dist_params Parameters for weight distribution.
#' @param coord_dist Distribution for coordinates.
#' @param w_dist Distribution for weights. 
#'
#' @return Data.frame of simulated data.
#' @export
#'
#' @examples
simulate_data <- function(n, w_dist_params = c(1, 200), coord_dist = "uniform", w_dist = "uniform") {
  
  lon <- switch(coord_dist, "uniform" = runif(n), "normal" = rnorm(n), "laplace" = rlaplace(n)) 
  lat <- switch(coord_dist, "uniform" = runif(n), "normal" = rnorm(n), "laplace" = rlaplace(n)) 
  w <- switch(w_dist,
              "uniform" = floor(runif(n, min = w_dist_params[1], max = w_dist_params[2] + 1)),
              "normal" = c(round(rnorm(n, mean = w_dist_params[1], sd = w_dist_params[2]))))
  w <- ifelse(w > 0, w, floor(runif(n, min = 1, max = quantile(w, probs = c(0.25)))))
  
  sim_data <- data.frame(id = 1:n,
                         lon = lon,
                         lat = lat,
                         max_load = w)
  return(sim_data)
}