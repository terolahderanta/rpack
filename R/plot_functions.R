#' Plot the clusters.
#' @param coords Coordinates of the data points.
#' @param weights The weigths of the objects in data.
#' @param clusters A vector of cluster assignments for each data point.
#' @param mu The parameters (locations) that define the k distributions.
#' @param title Set the title of the plot.
#' @param subtitle Set the subtitle of the plot.
#' @export
plot_clusters <- function(coords, weights, clusters, mu, title = "", subtitle = NULL){
  
  # x and y coordinates
  x <- coords[,1]
  y <- coords[,2]
  
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
