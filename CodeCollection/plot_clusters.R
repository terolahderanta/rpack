plot_clusters = function (coords, weights, clusters, centers, title = "", subtitle = NULL, outgroup_label=99, outgroup_legend = "outgroup") 
{
  
  coords = test_dat[,1:2]
  centers = as.matrix(centers)
  coords = as.matrix(coords)
  x = coords[, 1]
  y = coords[, 2]
  k = nrow(centers)
  w = test_dat$w
  centers = as.matrix(centers)
  help_clusters = clusters
  clusters = as.factor(clusters)
  cl_sizes = apply(X = t(1:k), MARGIN = 2, FUN = function(x) {
    sum(w[clusters == x])
  })
  cl_sizes = cl_sizes[help_clusters]
  cl_sizes[is.na(cl_sizes)] = sum(w[help_clusters == outgroup_label])
  cl_sizes[help_clusters == outgroup_label] = 
    paste(outgroup_legend, " (", cl_sizes[help_clusters == outgroup_label][1], ")", sep="")
  cl_sizes = factor(cl_sizes, 
                    levels = unique(cl_sizes)[order(nchar(unique(cl_sizes)), unique(cl_sizes))]) 
  plot = ggplot2::ggplot() + 
    ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = y, size = w, 
                                               colour = cl_sizes)) + 
    ggplot2::scale_size(range = c(2, 7), guide = FALSE) + 
    ggplot2::scale_color_manual(values = c_col[1:(k+1)], name = "Cluster sizes:") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) + 
    ggplot2::labs(x = "x", y = "y", title = NULL, subtitle = NULL) + 
    ggplot2::theme(legend.position = "right", axis.text.x = ggplot2::element_blank(), 
                   axis.text.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) + 
    ggplot2::geom_point(mapping = ggplot2::aes(x = centers[, 1], y = centers[, 2]), size = 3, 
                        col = "black", show.legend = FALSE, shape = 4, stroke = 3) + 
    ggplot2::ggtitle(subtitle) + 
    ggplot2::theme(
      plot.title = element_text(size=30, face="bold", margin = margin(10, 0, 10, 0)),
      axis.text.x = element_text(angle=0, size = 12),
      axis.text.y = element_text(angle = 0, size = 12),
      axis.title = element_text(size=14),
      legend.text = element_text(size=14), 
      legend.title = element_text(size = 16)
    ) +
    coord_fixed()
  
  return(plot)
  
}