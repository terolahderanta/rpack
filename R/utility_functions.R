c_col = c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
          "chocolate","darkgrey","darkred","yellow3","darkgreen","wheat3","magenta",
          "palegreen2","violetred","seagreen2","tomato4","steelblue1","royalblue",
          "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
          "coral3","beige","brown4","indianred1","lightgreen","orchid")

#' Calculate the mode of vector v
#' @param v A vector.
#' @return The mode.
#' @export
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Plot the clusters.
#' @param x x-coordinates.
#' @param y y-coordinates.
#' @param weights The weigths of the objects in data.
#' @param clusters A vector of cluster assignments for each data point.
#' @param mu The parameters (locations) that define the k distributions.
#' @param title Set the title of the plot.
#' @param subtitle Set the subtitle of the plot.
#' @export
plot_clusters <- function(x, y, weights, clusters, mu, title = "", subtitle = NULL){

  # The number of clusters
  k <- nrow(mu)

  # Changing clusters to factors
  clusters <- as.factor(clusters)

  # Cluster sizes
  cl_sizes <- apply(
    X = t(1:k),
    MARGIN = 2,
    FUN = function(x) { sum(weights[clusters == x]) }
  )

  # Plot the clusters with ggplot
  plot <-
    ggplot2::ggplot(data = NULL) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = x,
        y = y,
        size = weights,
        col = clusters
      )
    ) +
    ggplot2::scale_size(range = c(2, 7),        # Scale objects sizes
                        guide = FALSE) +
    ggplot2::scale_color_manual(  # Color theme for objects and legend title
      values = rep(c_col, times = 5),
      name = "Cluster sizes:",
      labels = cl_sizes
    ) +
    ggplot2::guides(                            # Point size in legend
      color = ggplot2::guide_legend(
        override.aes = list(size=5)
      )
    ) +
    ggplot2::labs(                          # Labels for axis and title
      x = "x1",
      y = "x2",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme(
      legend.position = "right",               # Legend position and removing ticks from axis
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = mu[,1], # Plot cluster centers
                             y = mu[,2]),
      size = 3,
      col = "black",
      show.legend = FALSE,
      shape = 4,
      stroke = 3
    )

  return(plot)
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

  cl <- stats::kmeans(X, X[C, ])
  cl$initial_centers <- X[C,]

  return(cl)
}
