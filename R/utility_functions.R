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
#' @param data FIXME
#' @param weights FIXME
#' @param clusters FIXME
#' @param mu FIXME
#' @param title Set the title of the plot.
#' @param subtitle Set the subtitle of the plot.
#' @export
plot_clusters <- function(data, weights, clusters, mu, title = "", subtitle = NULL){
  # Changing mu to data.frame
  mu <-  data.frame(mu)

  # The number of clusters
  k <- nrow(mu) #length(mu[, 1])

  # Changing clusters to factors
  clusters <- as.factor(clusters)

  # Cluster sizes
  cl_sizes <- apply(X = t(1:k), MARGIN = 2, FUN = function(x) {sum(weights[clusters == x])})

  # Plot the clusters with ggplot
  plot <-
    ggplot2::ggplot(data = NULL) +
    ggplot2::geom_point(mapping =
                          ggplot2::aes(x = data[, 1],
                                       y = data[, 2],
                                       size = weights,
                                       col = clusters)) +
    ggplot2::scale_size(range = c(2, 7), guide = FALSE) +    # Scale objects sizes
    ggplot2::scale_color_manual(values =                     # Color theme for objects and legend title
                                  rep(c_col,
                                      times = 5),
                                name = "Cluster sizes:",
                                labels = cl_sizes) +
    ggplot2::guides(color =                                 # Point size in legend
                      ggplot2::guide_legend(override.aes =
                                              list(size=5))) +
    ggplot2::labs(x = "x1", y = "x2",
                  title = title,
                  subtitle = subtitle) +      # Labels for axis and title
    ggplot2::theme(legend.position = "right",               # Legend position and removing ticks from axis
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = mu[, 1], # Plot cluster centers
                                               y = mu[, 2]),
                        size = 3,
                        col = "black",
                        show.legend = FALSE,
                        shape = 4,
                        stroke = 3)

  return(plot)
}

#' K-means++ algorithm FIXME: insert ref!
#' @param X FIXME
#' @param k FIXME
#' @export
kmpp <- function(X, k) {

  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)

  for (i in 2:k) {
    dm <- pracma::distmat(X, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }

  stats::kmeans(X, X[C, ])
}
