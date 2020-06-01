#' Capacitated location allocation  
#' 
#' Alternating algorithm for maximizing the optimization function in LA/clustering. Includes various constraints and extensions 
#' such as: Capacity limits, outliers, weights for the data points, different distance metrics, different
#' memberships, fixed points etc. 
#' 
#' @param coords A matrix or data.frame containing the coordinates.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param ranges Lower and upper limits for the clustering.
#' @param mob Mobility matrix.
#' @param capacity_weights Different weights for capacity limits.
#' @param d The distance function.
#' @param center_init Options to initialize center locations. Default is "random" and other choice is "kmpp". 
#' @param fixed_centers Predetermined center locations.
#' @param lambda Outlier-parameter.
#' @param place_to_point if TRUE, cluster centers will be placed to a point.
#' @param frac_memb If TRUE memberships are fractional.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param dist_mat Distance matrix for all the points. 
#' @param multip_centers Vector (n-length) defining how many centers a point is allocated to. 
#' @param parallel Should parallel computing be used?
#' @param print_output Print options, default is NULL. One option is "steps" which prints information about steps.
#' @return A list containing cluster allocation, cluster center and the current value of the objective function.
#' @keywords internal
capacitated_LA_mob <- function(coords,
                           weights,
                           k,
                           ranges,
                           mob,
                           lambda_mob = 1,
                           capacity_weights = weights,
                           d = euc_dist2, 
                           center_init = NULL,
                           fixed_centers = NULL, 
                           lambda = NULL, 
                           place_to_point = TRUE, 
                           frac_memb = FALSE, 
                           gurobi_params = NULL, 
                           dist_mat = NULL,
                           multip_centers = rep(1, nrow(coords)),
                           parallel = FALSE,
                           print_output = NULL,
                           lambda_fixed = NULL){
  
  # Number of objects in data
  n <- nrow(coords)
  
  # Number of fixed centers
  n_fixed <- ifelse(is.null(fixed_centers), 0, nrow(fixed_centers))
  
  # Cluster centers with fixed centers first
  if(n_fixed > 0){
    
    if(place_to_point){
      fixed_center_ids <- which(
        ((coords %>% pull(1)) %in% (fixed_centers %>% pull(1))) & 
          ((coords %>% pull(2)) %in% (fixed_centers %>% pull(2)))
      )
      
    } 
  } else {
    fixed_center_ids <- NULL
  }
  
  # Initialize centers
  switch(center_init,
         # Pick initial centers randomly
         random = {
           
           if(place_to_point){
             # Don't sample the points in fixed_centers
             sample_points <- which(!(1:n %in% fixed_center_ids))
             center_ids <- c(fixed_center_ids, sample(sample_points, k - n_fixed))
             centers <- coords[center_ids,]
           } else {
             center_ids <- NULL
             centers <- coords[sample(1:n, k),]
           }
         },
         
         # Pick initial centers with k-means++
         kmpp = {
           
           # Replicate data points according to their weight
           weighted_coords <- apply(coords, 2, function(x) rep(x, round(weights)))
           
           # Run k-means++
           init_kmpp <- kmpp(weighted_coords, k)
           
           if(place_to_point){
             
             # Select closest points to the kmpp-result as the centers
             center_ids <- apply(X = init_kmpp$centers, MARGIN = 1, FUN = function(x){
               temp_dist_kmpp <- 
                 apply(X = coords, MARGIN = 1, FUN = function(y){d(x,y)})
               which.min(temp_dist_kmpp)
             })
             
             # Do not choose the same point twice
             ifelse(duplicated(center_ids),
                    sample((1:n)[-center_ids],size = 1),
                    center_ids)
             
             # Select the centers according to ids
             centers <- coords[center_ids,]
             
           } else {
             centers <- init_kmpp$centers
             center_ids <- NULL
           }
         },
         
         stop("No such choice for center initialization! (rpack)")
  )
  
  # Maximum number of laps
  max_laps <- 100
  
  for (iter in 1:max_laps) {
    # Old mu is saved to check for convergence
    old_centers <- centers
    
    # Print detailed steps
    if(print_output == "steps"){
      cat("A-step... ")
      temp_time <- Sys.time()
    }
    
    # Clusters in equally weighted data (Allocation-step)
    temp_alloc <- allocation_step_mob(
      coords = coords,
      weights = weights,
      k = k,
      centers = centers,
      ranges = ranges,
      mob = mob,
      lambda_mob = lambda_mob,
      center_ids = center_ids,
      capacity_weights = capacity_weights,
      lambda = lambda,
      d = d,
      frac_memb = frac_memb,
      gurobi_params = gurobi_params,
      dist_mat = dist_mat,
      multip_centers = multip_centers
    )
    
    # Print detailed steps
    if(print_output == "steps"){
      cat(paste("Done! (", format(round(Sys.time() - temp_time, digits = 3), nsmall = 3), ")\n", sep = ""))
    }
    
    # Save the value of the objective function
    obj_min <- temp_alloc$obj_min
    
    # Print detailed steps
    if(print_output == "steps"){
      cat("L-step... ")
      temp_time <- Sys.time()
    }
    
    # Updating cluster centers (Parameter-step)
    temp_loc <- location_step_mob(
      coords = coords,
      weights = weights,
      k = k,
      assign_frac = temp_alloc$assign_frac,
      fixed_centers = fixed_centers,
      d = d,
      place_to_point = place_to_point,
      dist_mat = dist_mat,
      parallel = parallel,
      lambda_fixed = lambda_fixed
    )
    
    # Print detailed steps
    if(print_output == "steps"){
      cat(paste("Done! (", format(round(Sys.time() - temp_time, digits = 3), nsmall = 3), ")\n", sep = ""))
    }
    
    centers <- temp_loc$centers
    
    center_ids <- temp_loc$center_ids
    
    # Print a message showing that max number of iterations was reached
    if(iter == max_laps & print_output == "steps"){
      cat(paste("WARNING! Reached maximum number of LA-iterations! Returning the clustering from last lap...\n",sep = ""))
    }
    
    # If nothing is changing, stop
    if(all(old_centers == centers)) break
  }
  
  # Save the assignments
  assign_frac <- temp_alloc$assign_frac
  
  # Hard clusters from assign_frac
  clusters <- apply(assign_frac, 1, which.max)
  
  # Cluster 99 is the outgroup
  clusters <- ifelse(clusters == (k+1), 99, clusters)
  
  # Return cluster allocation, cluster center and the current value of the objective function
  return(list(clusters = clusters, centers = centers, obj = obj_min, assign_frac = assign_frac))
}

#' Update cluster allocations by minimizing the objective function.
#'
#' @param coords A matrix or data.frame containing the coordinates.
#' @param weights A vector of weights for each data point.
#' @param k The number of clusters.
#' @param centers The parameters (locations) that define the k distributions.
#' @param ranges Lower and upper limits for the clustering.
#' @param mob Mobility matrix.
#' @param center_ids Ids for the data points that are selected as centers.
#' @param capacity_weights Different weights for capacity limits.
#' @param lambda Outlier-parameter
#' @param d Distance function.
#' @param frac_memb If TRUE memberships are fractional.
#' @param gurobi_params A list of parameters for gurobi function e.g. time limit, number of threads.
#' @param dist_mat Distance matrix for all the points.
#' @param multip_centers Vector (n-length) defining how many centers a point is allocated to. 
#' @return New cluster allocations for each object in data and the maximum of the objective function.
#' @keywords internal
allocation_step_mob <- function(coords,
                            weights,
                            k,
                            centers,
                            ranges,
                            mob,
                            lambda_mob = 1,
                            center_ids = NULL,
                            capacity_weights = weights,
                            lambda = NULL,
                            d = euc_dist2,
                            frac_memb = FALSE,
                            gurobi_params = NULL,
                            dist_mat = NULL,
                            multip_centers = rep(1, nrow(coords))) {
  
  # Number of objects in data
  n <- nrow(coords)
  
  # Is there an outgroup cluster
  is_outgroup <- !is.null(lambda)
  
  # Number of range groups
  if(is.vector(ranges)){
    ranges <- matrix(data = ranges, nrow = 1, ncol = 2)
    g <- 1
  } else {
    g <- nrow(ranges) 
  }
  
  # Number of decision variables
  n_decision <- n * k + 
    # More than one range groups
    ifelse(g > 1, k * g, 0) +
    # Outliers
    ifelse(is_outgroup, n, 0)
  
  # Calculate the distances to centers (matrix C)
  if(is.null(dist_mat) | length(center_ids) == 0){
    
    C <- matrix(0, ncol = k, nrow = n)
    for(i in 1:k){
      C[,i] <- apply(coords, MARGIN = 1, FUN = d, x2 = centers[i,])
    }
    
  } else {
    # Read distances from distance matrix
    C <- dist_mat[,center_ids]
  }
  
  # Use weighted distances
  C <- C * weights
  
  # Gurobi model
  model <- list()
  
  if(g == 1){
    
    # Constraints for the upper and lower limit
    const_LU <- rbind(
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = rep(1:k, each = n),
        j = 1:(n * k),
        x = rep(capacity_weights, k)
      ),
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = rep(1:k, each = n),
        j = 1:(n * k),
        x = rep(capacity_weights, k)
      )
    )
    
    # Add the constraints to the model
    model$A <- rbind(Matrix::spMatrix(
      nrow = n,
      ncol = n_decision,
      i = rep(1:n, times = ifelse(is_outgroup, k + 1, k)),
      j = rep(1:n_decision),
      x = rep(1, n_decision)
    ),
    const_LU)
    
    # Right hand side values
    model$rhs <- c(multip_centers,
                   rep(ranges[1, 2], k),
                   rep(ranges[1, 1], k))
    
    # Model sense
    model$sense      <- c(rep('=', n), 
                          rep('<', k), 
                          rep('>', k))
    
    
  } else {
    # Large number
    M <- 1000
    
    # In constraint matrices: first rangegroups, then clusters
    
    # Constraints for the first lower and upper limit
    const_LU <- rbind(
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = c(rep(1:k, each = n), 1:k),
        j = c(1:(n * k), (n * k) + 1:k),
        x = c(rep(capacity_weights, k), rep(M, k))
      ),
      Matrix::spMatrix(
        nrow = k,
        ncol = n_decision,
        i = c(rep(1:k, each = n), 1:k),
        j = c(1:(n * k), (n * k) + 1:k),
        x = c(rep(capacity_weights, k), rep(-M, k))
      )
    )
    
    # Right hand side values for the first capacity constraints
    rhs_LU <- c(rep(ranges[1, 2] + M, k),
                rep(ranges[1, 1] - M, k))
    
    # Model sense for the first capacity constraints
    sense_LU     <- c(rep('<', k),
                      rep('>', k))
    
    # Constraints, rhs and sense for the rest of the lower and upper limits
    for(i in 2:g){
      const_LU <- rbind(
        const_LU,
        Matrix::spMatrix(
          nrow = k,
          ncol = n_decision,
          i = c(rep(1:k, each = n), 1:k),
          j = c(1:(n * k), (n * k + k * (i - 1)) + 1:k),
          x = c(rep(capacity_weights, k), rep(M, k))
        ),
        Matrix::spMatrix(
          nrow = k,
          ncol = n_decision,
          i = c(rep(1:k, each = n), 1:k),
          j = c(1:(n * k), (n * k + k * (i - 1)) + 1:k),
          x = c(rep(capacity_weights, k), rep(-M, k))
        )
      )
      
      rhs_LU <- c(rhs_LU,
                  rep(ranges[i, 2] + M, k),
                  rep(ranges[i, 1] - M, k))
      
      sense_LU <- c(sense_LU,
                    rep('<', k),
                    rep('>', k))
      
    }
    
    # Constraints for the cluster size group
    const_group <- Matrix::spMatrix(
      nrow = k,
      ncol = n_decision,
      i = rep(1:k, each = g),
      j = c((n * k) + 1:(k * g)),
      x = rep(1, k * g)
    )
    
    # Add all constraints to the model
    model$A <- rbind(
      Matrix::spMatrix(
        nrow = n,
        ncol = n_decision,
        i = rep(1:n, times = ifelse(is_outgroup, k + 1, k)),
        j = rep(1:(n * k)),
        x = rep(1, n * k)
      ),
      const_LU,
      const_group
    )
    
    # Right hand side values (multiple membership, upper and lower limits, cluster groups)
    model$rhs <- c(multip_centers,
                   rhs_LU,
                   rep(1, k))
    
    # Model sense
    model$sense <- c(rep('=', n),
                     sense_LU,
                     rep('=', k))
  }
  
  # Indeces of nonzero values in the mobility matrix
  mob_nonzero_ind <- which(mob != 0, arr.ind = TRUE)
  
  # Nonzero values of the mobility matrix
  mob_nonzero_val <- mob[mob_nonzero_ind]
  
  # Transform the mobility matrix into a sparse matrix given to gurobi
  mob_Q <- Matrix::sparseMatrix(i = rep(mob_nonzero_ind[, 1], k) + n * rep(0:(k-1), each = length(mob_nonzero_ind[, 1])), 
                        j = rep(mob_nonzero_ind[, 2], k) + n * rep(0:(k-1), each = length(mob_nonzero_ind[, 2])),
                        x = rep(mob_nonzero_val, k),dims = c(n_decision, n_decision))
  
  # Quadratic objective matrix
  Q <- matrix(0.5, ncol = n_decision, nrow = n_decision) + diag(0.5, ncol = n_decision, nrow = n_decision)
  
  model$Q <- -lambda_mob * Q * mob_Q
  
  
  # Objective function
  obj_fn <- c(c(C),
              switch(g > 1, rep(0, k * g), NULL),
              switch(is_outgroup, lambda * weights, NULL))
  
  model$obj <- obj_fn
  
  # Minimization task
  model$modelsense <- 'min'
  
  # B = Binary, C = Continuous
  model$vtype <- ifelse(frac_memb, 'C', 'B')
  
  # Using timelimit-parameter to stop the optimization if time exceeds 10 minutes
  # and disabling the print output from gurobi.
  if(is.null(gurobi_params)){
    gurobi_params <- list()
    gurobi_params$TimeLimit <- 600
    gurobi_params$OutputFlag <- 0  
  }
  
  # Solving the linear program
  result <- gurobi::gurobi(model, params = gurobi_params)
  
  # Send error message if the model was infeasible
  if(result$status == "INFEASIBLE") {stop("Model was infeasible! (rpack)")}
  
  # Returns the assignments
  assign_frac <- Matrix::Matrix(matrix((result$x)[1:(ifelse(is_outgroup, n * k + n, n * k))],
                                       ncol = ifelse(is_outgroup, k + 1, k)), sparse = TRUE)
  
  # Returns the value of the objective function
  obj_min <- round(result$objval, digits = 5)
  
  # Clear space
  rm(model, result)
  
  return(list(assign_frac = assign_frac,
              obj_min = obj_min))
}


#' Updates the parameters (centers) for each cluster.
#'
#' @param coords A matrix or data.frame containing the coordinates.
#' @param weights The weights of the data points.
#' @param k The number of clusters.
#' @param assign_frac A vector of cluster assignments for each data point.
#' @param fixed_centers Predetermined center locations.
#' @param d The distance function.
#' @param place_to_point if TRUE, cluster centers will be placed to a point.
#' @param dist_mat Distance matrix for all the points.
#' @param parallel Logical indicator to use parallel computing.
#' @return New cluster centers.
#' @keywords internal
location_step_mob <- function(coords,
                          weights,
                          k,
                          assign_frac,
                          fixed_centers = NULL,
                          d = euc_dist2,
                          place_to_point = TRUE,
                          dist_mat = NULL,
                          parallel = FALSE,
                          lambda_fixed = NULL) {
  
  # Number of fixed centers
  n_fixed <- ifelse(is.null(fixed_centers), 0, nrow(fixed_centers))
  
  # Initialization of cluster center matrix
  centers <- matrix(0, nrow = k, ncol = ncol(coords))
  
  # Add fixed centers first
  if(n_fixed > 0){
    centers[1:n_fixed,] <- fixed_centers 
  }
  
  if(place_to_point){
    # Initialization of cluster id vector
    center_ids <- rep(0, k)
    
    # Add fixed centers first
    if(n_fixed > 0){
      center_ids <-c(which(
        ((coords %>% pull(1)) %in% (fixed_centers %>% pull(1))) & 
          ((coords %>% pull(2)) %in% (fixed_centers %>% pull(2)))
      ), rep(0, k - n_fixed))
    }
  }
  
  # Use parallel computing
  if (parallel) {
    # Setup parallel backend to use all but one processor
    cores <- detectCores()
    cl <- makeCluster(cores[1] - 1)
    registerDoParallel(cl)
    
    if (n_fixed > 0) {
      stop("Can't use fixed centers and parallel together.... yet (rpack)")
    }
    
    # Update center of each cluster
    if (place_to_point) {
      
      center_ids <- rep(0,k)
      
      center_ids <-
        foreach(i = (n_fixed + 1):k, .combine = rbind, .packages = "Matrix") %dopar% {
          # Compute medoids only with points that are relevant in the cluster i
          relevant_cl <- assign_frac[, i] > 0.001
          
          # Computing medoid ids for cluster i
          temp_center_id <- c(medoid_dist_mat(dist_mat = dist_mat,
                                              ids = which(relevant_cl),
                                              w = weights))
          
          # rbind the temp_center
          temp_center_id
        }
      
      
      # Decide centers from the ids
      centers <- coords[center_ids, ]
      
    } else {
      centers <-
        foreach(i = (n_fixed + 1):k, .combine = rbind, .packages = "Matrix") %dopar% {
          # Check whether euc_dist or euc_dist2 is used
          if (d(0, 2) == 2) {
            # Weighted median
            temp_center <-
              as.matrix(Gmedian::Weiszfeld(coords, weights = weights * assign_frac[, i])$median)
            
          } else if (d(0, 2) == 4) {
            
            # Weighted mean
            temp_center <-
              colSums(coords * weights * assign_frac[, i]) / sum(assign_frac[, i] * weights)
          }
          
          # rbind the temp_center
          matrix(temp_center, ncol = 2)
        }
      
      center_ids <- NULL
    }
    
    # Stop cluster
    stopCluster(cl)
    
  } else {
    # Update center for each cluster
    if (place_to_point) {
      
      # If fixed servers are allowed to be released
      if(!is.null(lambda_fixed)) {
        
        # Potential new centers
        potential_center_ids <- rep(0,n_fixed)
        
        # Ids for the fixed servers
        fixed_center_ids <- which(
          ((coords %>% pull(1)) %in% (fixed_centers %>% pull(1))) & 
            ((coords %>% pull(2)) %in% (fixed_centers %>% pull(2)))
        )
        
        # Most optimal center location for the fixed center clusters
        for(i in 1:n_fixed){
          # Compute medoids only with points that are relevant in the cluster i
          relevant_cl <- assign_frac[, i] > 0.001
          
          relevant_ids <- which(relevant_cl)
          
          # Computing medoid ids for cluster i
          potential_center_ids[i] <- medoid_dist_mat(dist_mat = dist_mat,
                                                     ids = relevant_ids,
                                                     w = weights)
          
          # Combined distance to the potential center
          wdist_pot_center <- sum(
            dist_mat[relevant_ids, potential_center_ids[i]] * 
              weights[relevant_ids]
          )
          
          # Combined distance to the fixed center
          wdist_fixed_center <- sum(
            dist_mat[relevant_ids, fixed_center_ids[i]] * 
              weights[relevant_ids]
          )
          
          # TODO: add lambda_fixed to here
          if(wdist_fixed_center > wdist_pot_center + lambda_fixed){
            center_ids[i] <- potential_center_ids[i]
          }
          
        }
        
        # Decide the rest of the center locations
        for (i in (n_fixed + 1):k) {
          
          # Compute medoids only with points that are relevant in the cluster i
          relevant_cl <- assign_frac[, i] > 0.001
          
          # Computing medoid ids for cluster i
          center_ids[i] <- medoid_dist_mat(dist_mat = dist_mat,
                                           ids = which(relevant_cl),
                                           w = weights)
        }
        
        # Decide centers from the ids
        centers <- coords[center_ids,]
        
      } else {
        
        for (i in (n_fixed + 1):k) {
          
          # Compute medoids only with points that are relevant in the cluster i
          relevant_cl <- assign_frac[, i] > 0.001
          
          # Computing medoid ids for cluster i
          center_ids[i] <- medoid_dist_mat(dist_mat = dist_mat,
                                           ids = which(relevant_cl),
                                           w = weights)
        }
        
        # Decide centers from the ids
        centers <- coords[center_ids,]
      }
    } else {
      for (i in (ifelse(n_fixed > 0, n_fixed + 1, 1)):k) {
        # Check whether euc_dist or euc_dist2 is used
        if (d(0, 2) == 2) {
          
          w_assign <- weights * assign_frac[, i]
          
          # If only one point belongs to the cluster, the Weiszfeld algorithm won't work
          if(sum(w_assign > 0) == 1){
            centers[i, ] <- coords %>% 
              slice(which(w_assign > 0)) %>% 
              unlist(., use.names=FALSE)
            
          } else {
            # Weighted median
            weiszfeld <- Gmedian::Weiszfeld(coords, weights = w_assign)$median 
            
            centers[i, ] <- weiszfeld
            
          }
          
        } else if (d(0, 2) == 4) {
          # Weighted mean
          centers[i, ] <-
            colSums(coords * weights * assign_frac[, i]) / sum(assign_frac[, i] * weights)
        }
      }
      
      center_ids <- NULL
    }
    
  }
  return(list(centers = centers, center_ids = center_ids))
}
