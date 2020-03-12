#' Different colors for the plots
#' @export
c_col <- c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
          "chocolate","darkred","yellow3","darkgreen","bisque4","magenta",
          "royalblue","tomato4","steelblue1",
          "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
          "coral3","beige","brown4","indianred1","lightgreen","orchid")

#' Calculate the mode of vector v
#' @param v A vector.
#' @return The mode.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Calculate the Euclidean distance between two points
#' @param x1 1. point.
#' @param x2 2. point.
#' @return Euclidean distance.
#' @export
euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#' Calculate the squared Euclidean distance between two points
#' @param x1 1. point.
#' @param x2 2. point.
#' @return Squared Euclidean distance.
#' @export
euc_dist2 <- function(x1, x2) sum((x1 - x2) ^ 2)

#' Calculate the medoid of the data points
#' @param data A data.frame.
#' @param w Weights of the data points.
#' @param d A distance metric.
#' @export
#' @return The medoid.
medoid <- function(data,
                   w = rep(1, nrow(data)),
                   d = euc_dist2) {
  n <- nrow(data)
  if (n < 1) {
    stop("Tried to calculate medoid from zero number of points! (rpack)")
  }
  if (n == 1) {
    return(data[1, ])
  }
  
  w_dists <- sapply(
    1:n,
    FUN = function(x) {
      sum(w[-x] * apply(data[-x, ], 1, FUN = d, x2 = data[x, ]))
    }
  )
  
  return(data[which.min(w_dists), ])
}

#' Calculate the medoid from distance matrix
#' @param dist_mat Distance matrix for the data points.
#' @param ids Ids for the points in distance matrix. Uses all of the points by default.  
#' @param w Weights of the data points.
#' @export
#' @return The id for the medoid.
medoid_dist_mat <- function(dist_mat,
                            ids = 1:nrow(dist_mat),
                            w = rep(1, nrow(dist_mat))) {
  
  # Exceptions
  n <- nrow(dist_mat)
  if (n < 1 | length(ids) == 0) {
    stop("Tried to calculate medoid from zero number of points! (rpack)")
  }
  if (n == 1 | length(ids) == 1) {
    return(ids[1])
  }
  
  # Weighted distances from the given set of points
  wdists <- dist_mat[ids, ids] * w[ids]
  
  # Calculate column sums
  wdist_to_centers <- colSums(wdists)
  
  return(ids[which.min(wdist_to_centers)])
}

#' Kmeans++
#'
#' Implementation of the K-means++ algorithm. Whereas normal kmeans selects all the initial center
#' cluster centers randomly, kmeans++ randomly selects only the first center. For each
#' consecutive center, the probability of selection is weighed by the distance to already selected
#' centers.
#'
#' Implementation adapted from one by Hans Werner (https://stat.ethz.ch/pipermail/r-help/2012-January/300051.html).
#'
#' See following article for more information on kmeans++:
#'
#'   Arthur, D., and Vassilvitskii, S. "k-means++: The advantages of careful seeding."
#'   Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms. Society
#'   for Industrial and Applied Mathematics, 2007.
#'
#' @param X A matrix or a data frame containing the objects, one per row.
#' @param k Number of clusters.
#' @export
kmpp <- function(X, k) {

  if (!is.matrix(X)) X <- as.matrix(X)  # X must be a matrix

  n <- nrow(X)
  ncoords <- ncol(X)
  C <- numeric(k)                # initialize centers to zero
  C[1] <- sample(1:n, size = 1)  # select first element randomly

  for (i in 2:k) {
    dm <- pracma::distmat(X, matrix(X[C, ], ncol = ncoords))
    pr <- apply(dm, 1, min)
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  if(length(unique(C)) == k){
    cl <- stats::kmeans(X, X[C, ])  
    cl$initial_centers <- X[C,]
  } else {
    cl <- kmpp(X, k)
  }
  #cl <- stats::kmeans(X, X[C, ])
  #cl$initial_centers <- X[C,]

  return(cl)
}
