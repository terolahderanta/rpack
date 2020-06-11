#' Plot the clusters with a convex hull.
#'
#' @param coords Coordinates of the data points.
#' @param weights The weights of the objects in data.
#' @param clusters A vector of cluster assignments for each data point.
#' @param centers The centers of the clusters.
#' @param assing_frac Assignments to the clusters in case of a fractional membership.
#' @param fixed Coordinates of fixed centers.
#' @param title Title for the plot.
#' @param alpha Level of transparency.
#' @param frac_memb Is the clustering done with fractional membership. 
#'
#' @return ggplot object.
#' @export
plot_hull <- function(coords, weights, clusters, centers, assign_frac = NULL, fixed = NULL, title = "", alpha = 0.4, frac_memb = FALSE) {
  
  font_size <- 15
  
  # The number of points
  n <- nrow(coords)
  
  # The number of clusters
  k <- nrow(centers)
  
  # Transform centers into a matrix
  centers <- as.matrix(centers)
  
  # Find the convex hull of the points being plotted
  hull <- list()
  
  if(!frac_memb) {
    hull <- cbind(coords[clusters == 1, ] %>% slice(chull(x, y)), cl = 1)
    for (i in 2:k) {
      hull <- rbind(hull, cbind(coords[clusters == i, ] %>% slice(chull(x, y)), cl = i))
    }
  } else {
    
    if(is.null(assign_frac)){stop("Provide assignment matrix when using fractional membership! (rpack)")}
    
    secondary <- rep(0, n)
    for (i in 1:n) {
      # Find the second biggest assignment for point i
      sec_max <- sort(assign_frac[i,], partial=k-1)[k-1]
      
      # If second biggest assignment is bigger than 0.001, then it is secondary assignment
      if(sec_max > 0.01){
        secondary[i] <- which(assign_frac[i,] == sec_max)
      }
    }
    hull <- cbind(coords[clusters == 1 | secondary == 1, ] %>% slice(chull(x, y)), cl = 1)
    for (i in 2:k) {
      hull <- rbind(hull, cbind(coords[clusters == i | secondary == i, ] %>% slice(chull(x, y)), cl = i))
    }
    
  }
  
  mu <- data.frame(x = centers[,1],
                   y = centers[,2])
  if(!is.null(fixed)) {
    fixed_mu <- data.frame(x = fixed[,1],
                           y = fixed[,2])
  }
  
  # Calculate cluster sizes
  if(frac_memb){
    cl_sizes <- apply(
      X = t(1:k),
      MARGIN = 2,
      FUN = function(x) { sum(weights*assign_frac[,x]) }
    )
  } else {
    cl_sizes <- apply(
      X = t(1:k),
      MARGIN = 2,
      FUN = function(x) { sum(weights[clusters == x]) }
    )
  }
  
  # Order the clusters from smallest to larger
  order_incr <- order(cl_sizes)
  
  ggplot() + 
    geom_point(data = NULL, 
               aes(x = coords %>% pull(1), 
                   y = coords %>% pull(2),
                   size = weights,
                   color = factor(clusters, levels = as.factor((1:k)[order_incr])))) + 
    scale_size(range = c(1, 3),
               guide = FALSE) +
    #scale_x_continuous(limits = c(25.39, 25.61), expand = c(0, 0)) +
    #scale_y_continuous(limits = c(64.95, 65.09), expand = c(0, 0)) +
    xlab("") + 
    ylab("") +
    ggtitle(title) + 
    scale_color_manual(values = rep(c_col, times = 5),
                       name = "Cluster sizes:",
                       labels = cl_sizes[order_incr]) + 
    #scale_color_discrete(name = "Cluster sizes:",
    #                     labels = cl_sizes) +
    guides(color = ggplot2::guide_legend(override.aes = list(size=4))) +
    geom_polygon(data = hull, aes(x = x, 
                                  y = y,
                                  group = cl,
                                  fill = factor(cl, levels = as.factor((1:k)[order_incr]))),
                 alpha = alpha) +
    scale_fill_manual(guide = FALSE, values = rep(c_col, 5)) +
    geom_point(data = mu,
               mapping = aes(x = x,
                             y = y),
               size = 2,
               col = "black",
               show.legend = FALSE,
               shape = 4,
               stroke = 1
    ) + 
    theme(text = element_text(size=font_size)) +
    switch(is.null(fixed) + 1, geom_point(data = fixed_mu,
                                          mapping = aes(x = x,
                                                        y = y),
                                          size = 2,
                                          col = "black",
                                          show.legend = FALSE,
                                          shape = 3,
                                          stroke = 1
    ), NULL)
}

#' Plot the clusters.
#' @param coords Coordinates of the data points.
#' @param weights The weights of the objects in data.
#' @param clusters A vector of cluster assignments for each data point.
#' @param centers The centers of the clusters.
#' @param frac_memb Is the clustering done with fractional membership. 
#' @param title Set the title of the plot.
#' @param subtitle Set the subtitle of the plot.
#' @export
plot_clusters <- function(coords, weights, clusters, centers, frac_memb = FALSE, title = "", subtitle = NULL){
  
  # x and y coordinates
  coords <- as.matrix(coords)
  x <- coords[,1]
  y <- coords[,2]
         
  # The number of clusters
  k <- nrow(centers)
  
  # Transform centers into a matrix
  centers <- as.matrix(centers)
  
  # Changing clusters to factors
  clusters <- as.factor(clusters)
  
  # Calculate cluster sizes
  if(frac_memb){
    cl_sizes <- apply(
      X = t(1:k),
      MARGIN = 2,
      FUN = function(x) { sum(weights*assign_frac[,x]) }
    )
  } else {
    cl_sizes <- apply(
      X = t(1:k),
      MARGIN = 2,
      FUN = function(x) { sum(weights[clusters == x]) }
    )
  }
  
  # Order the clusters from smallest to larger
  order_incr <- order(cl_sizes)
  
  # Plot the clusters with ggplot
  plot <-
    ggplot2::ggplot(data = NULL) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = x,
        y = y,
        size = weights,
        col = factor(clusters, levels = as.factor((1:k)[order_incr]))
        )
    ) +
    ggplot2::scale_size(range = c(2, 7),        # Scale objects sizes
                        guide = FALSE) +
    ggplot2::scale_color_manual(  # Color theme for objects and legend title
      values = rep(c_col, times = 5),
      name = "Cluster sizes:",
      labels = cl_sizes[order_incr]
    ) +
    ggplot2::guides(                            # Point size in legend
      color = ggplot2::guide_legend(
        override.aes = list(size=5)
      )
    ) +
    ggplot2::labs(                          # Labels for axis and title
      x = "x",
      y = "y",
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
      mapping = ggplot2::aes(x = centers[,1], # Plot cluster centers
                             y = centers[,2]),
      size = 3,
      col = "black",
      show.legend = FALSE,
      shape = 4,
      stroke = 3
    )
  
  return(plot)
}
