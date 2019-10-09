#' Probabilistic Clustering
#'
#' Alternating algorithm for maximizing the joint density.
#'
#' @param data A matrix or data.frame whose rows contain the objects to be clustered.
#' @param weights A vector of weights for each data point. The weights will be rounded to nearest integers.
#' @param k The number of clusters.
#' @param capacity_weights Different weights for capacity limits.
#' @param d The distance function.
#' @param init_mu A matrix whose rows indicate the initial cluster centers.
#' @param prior_dist Prior distribution for cluster sizes. Possible values are "uniform" and "normal".
#' @param range Range of the uniform prior distribution as a vector c(lower, upper).
#' @param sigma Standard deviation of the normal prior.
#' @param lambda Outgroup-parameter.
#' @param divide_objects If TRUE, objects can be divided to multiple clusters
#' @param use_gurobi If TRUE, gurobi solver will be used in the optimization task
#' @param place_to_point if TRUE, cluster centers will be placed to a point.
#' @param fixed_mu Predetermined center locations.
#' @param frac_memb If TRUE memberships are fractional.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @return A list containting the new cluster allocations for each object in data,
#' the new cluster center locations and maximum of the objective function.
#' @export
prob_clust <- function(data,
                       weights,
                       k,
                       capacity_weights = weights,
                       d = euc_dist2,
                       init_mu = NULL,
                       prior_dist = "uniform",
                       range = NULL,
                       sigma = NULL,
                       lambda = NULL,
                       divide_objects = FALSE,
                       use_gurobi = TRUE,
                       place_to_point = TRUE,
                       fixed_mu = NULL,
                       frac_memb = FALSE,
                       gurobi_params = NULL) {
  
  # Check arguments
  assertthat::assert_that(is.matrix(data) || is.data.frame(data), msg = "data must be a matrix or a data.frame!")
  if (is.matrix(data)) data <- tibble::as_tibble(data)                   # convert to tibble for consistency
  
  assertthat::assert_that(nrow(data) >= k, msg = "must have at least k data points!")
  assertthat::assert_that(is.numeric(weights), msg = "weight must be an numeric vector!")
  assertthat::assert_that(length(weights) == nrow(data), msg = "data and weight must have the same number of rows!")
  assertthat::assert_that(is.numeric(k), msg = "k must be a numeric scalar!")
  assertthat::assert_that(length(k) == 1, msg = "k must be a numeric scalar!")
  
  if(!purrr::is_null(init_mu)) assertthat::assert_that(is.matrix(init_mu))
  if(!purrr::is_null(range)) {
    assertthat::assert_that(is.numeric(range))
    assertthat::assert_that(length(range) == 2)
  }
  if(!purrr::is_null(sigma)) assertthat::is.number(sigma)
  if(!purrr::is_null(lambda)) assertthat::is.number(lambda)
  
  assertthat::assert_that(is.logical(divide_objects), msg = "divide_objects must be TRUE or FALSE!")
  assertthat::assert_that(is.logical(use_gurobi), msg = "use_gurobi must be TRUE or FALSE!")
  
  # Create initial values for mu, if init_mu is not defined
  if(is.null(init_mu)){
    weighted_data <- apply(data, 2, function(x) rep(x, weights))  # replicate data points according to their weight
    init_kmpp <- kmpp(weighted_data, k)
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
      # Lower and upper limit given by user
      L <- range[1]
      U <- range[2]
    }
    
    if(use_gurobi){
      output_list <-
        prob_clust_gurobi(
          data = data,
          weights = weights,
          k = k,
          init_mu = init_mu,
          L = L,
          U = U,
          capacity_weights = capacity_weights,
          d = d,
          fixed_mu = fixed_mu,
          lambda = lambda,
          place_to_point = place_to_point,
          frac_memb = frac_memb,
          gurobi_params = gurobi_params
        )  
    } else {
      # Call function prob_clust_simple
      output_list <-
        prob_clust_uniform(
          data = data,
          weights = weights,
          k = k,
          init_mu = init_mu,
          L = L, U = U,
          d = d,
          lambda = lambda
        )
    }
  } else if(prior_dist == "normal"){
    
    # Initializing sigma
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
    prob <- stats::dnorm(cl_size, mean = pr_mean, sd = sigma)
    
    # Call function prob_clust_prior
    output_list <-
      prob_clust_prior(
        data = data,
        weights = weights,
        k = k,
        init_mu = init_mu,
        prior_cl_sizes = cl_size,
        prior_prob = prob,
        divide_objects = divide_objects
      )
    
  } else {
    stop("prior_dist must be 'uniform' or 'normal'.")
  }
  
  return(output_list)
}
