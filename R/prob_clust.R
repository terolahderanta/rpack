#' Probabilistic clustering algorithm.
#'
#' @param data A matrix or data.frame containing the data.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param init_mu Parameters (locations) that define the k distributions.
#' @param prior_dist Prior distribution for cluster sizes. Possible values are "uniform" and "normal".
#' @param range Range of the uniform prior distribution as a vector c(lower, upper).
#' @param sigma Standard deviation of the normal prior.
#' @param divide_objects If TRUE, objects can be divided to multiple clusters
#' @return A list containting the new cluster allocations for each object in data, the new cluster center locations and maximum of the objective function.
prob_clust <- function(data, weights, k, init_mu = NULL, prior_dist = "uniform", range = NULL, sigma = NULL, lambda = NULL, divide_objects = FALSE){
  # Creates initial values for mu, if init_mu is not defined
  if(is.null(init_mu)){
    init_kmpp <- kmpp(cbind(rep(data[, 1], weights), rep(data[, 2], weights)), k)
    init_mu <- init_kmpp$centers
  }

  # In case of uniform prior
  if(prior_dist == "uniform"){
    if(is.null(range)){
      # Mean for prior
      pr_mean <- round(sum(weights) / k)

      # Width of uniform prior
      pr_width <- round(max(weights) * 2)

      # Lower und upper limit for cluster size
      L <- (pr_mean - pr_width)
      U <- (pr_mean + pr_width)
    } else {
      L <- range[1]
      U <- range[2]
    }

    # Call function prob_clust_simple
    output_list <- prob_clust_simple(data = data, weights = weights, k = k, init_mu = init_mu, L = L, U = U, lambda = lambda)

  } else if(prior_dist == "normal"){

    if(is.null(sigma)){
      sigma <- max(weights) / 2
    }

    # Mean of the normal
    pr_mean <- round(sum(weights) / k)

    # Max width for prior
    pr_width <- round(sigma * 4)

    # All possible cluster sizes
    cl_size <- (pr_mean - pr_width):(pr_mean + pr_width)

    # Normal prior probabilities
    prob <- dnorm(cl_size, mean = pr_mean, sd = sigma)

    # Call function prob_clust_prior
    output_list <- prob_clust_prior(data = data, weights = weights, k = k, init_mu = init_mu,
                                     prior_cl_sizes = cl_size, prior_prob = prob, divide_objects = divide_objects)

  } else {
    stop("Prior distribution must be uniform or normal")
  }

  return(output_list)
}
