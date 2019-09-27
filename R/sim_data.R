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
simulate_data <- function(n, w_dist_params = c(1, 100), coord_dist = "uniform", w_dist = "uniform") {
  
  x <- switch(coord_dist, "uniform" = runif(n), "normal" = rnorm(n), "laplace" = rlaplace(n)) 
  y <- switch(coord_dist, "uniform" = runif(n), "normal" = rnorm(n), "laplace" = rlaplace(n)) 
  w <- switch(w_dist,
              "uniform" = floor(runif(n, min = w_dist_params[1], max = w_dist_params[2] + 1)),
              "normal" = c(round(rnorm(n, mean = w_dist_params[1], sd = w_dist_params[2]))))
  w <- ifelse(w > 0, w, floor(runif(n, min = 1, max = quantile(w, probs = c(0.25)))))
  
  sim_data <- data.frame(id = 1:n,
                         x = x,
                         y = y,
                         w = w)
  return(sim_data)
}

random_sigma <- function(){
  A <- matrix(runif(4, min = 0.1, max = 0.8),ncol = 2, nrow = 2)
  Sigma <- t(A) %*% A 
  Sigma[1,1] <- Sigma[1,1] + runif(1)
  Sigma[2,2] <- Sigma[2,2] + runif(1)
  return(Sigma)
}

#' Simulate from mixed normal distribution.
#'
#' @param n Total number of points simulated.
#' @param k Number of different normals.
#' @param w_dist_params Parameters for weight distribution
#' @param w_dist Distribution for weights.
#'
#' @return Data.frame of simulated data.
#' @export
#'
#' @examples
simulate_normal_mixed <- function(n, k, w_dist_params = c(1, 100), w_dist = "uniform") {
  
  n_sub <- round(n/k)
  
  # Mus for the normal distributions
  mu <- matrix(runif(2 * k, min = -10, max = 10), ncol = 2, nrow = k)
  
  # Sigmas for the normal distributions
  coords <- mvrnorm(n = n_sub,
                    mu = mu[1,],
                    Sigma = random_sigma())
  orig_group <- rep(1, n_sub)
  for (i in 2:k) {
    coords <- rbind(coords, mvrnorm(n = n_sub,
                                    mu = mu[i,],
                                    Sigma =  random_sigma()))
    orig_group <- c(orig_group, rep(i, n_sub))
  }
  
  # Weights
  w <- floor(runif(n, min = w_dist_params[1], max = w_dist_params[2] + 1))
  
  return(data.frame(x = coords[,1],
                    y = coords[,2],
                    w = w,
                    orig_group = as.factor(orig_group)))
}


#' Simulate from mixed normal distribution.
#'
#' @param n Total number of points simulated.
#' @param k Number of different normals.
#' @param w_dist_params Parameters for weight distribution
#' @param w_dist Distribution for weights.
#'
#' @return Data.frame of simulated data.
#' @export
#'
#' @examples
simulate_laplace_mixed <- function(n, k, w_dist_params = c(1, 100), w_dist = "uniform") {
  
  n_sub <- round(n/k)
  
  # Mus for the normal distributions
  mu <- matrix(runif(2 * k, min = -10, max = 10), ncol = 2, nrow = k)
  
  # Sigmas for the normal distributions
  coords <- rmvl(n = n_sub,
                    mu = mu[1,],
                    Sigma = random_sigma())
  orig_group <- rep(1, n_sub)
  for (i in 2:k) {
    coords <- rbind(coords, rmvl(n = n_sub,
                                    mu = mu[i,],
                                    Sigma =  random_sigma()))
    orig_group <- c(orig_group, rep(i, n_sub))
  }
  
  # Weights
  w <- floor(runif(n, min = w_dist_params[1], max = w_dist_params[2] + 1))
  
  return(data.frame(x = coords[,1],
                    y = coords[,2],
                    w = w,
                    orig_group = as.factor(orig_group)))
}
