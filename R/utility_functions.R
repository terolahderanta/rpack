c_col = c("blue","red","green","orange","hotpink","cyan","yellowgreen","purple",
          "chocolate","darkgrey","darkred","yellow3","darkgreen","wheat3","magenta",
          "palegreen2","violetred","seagreen2","tomato4","steelblue1","royalblue",
          "seagreen4","orangered","darkblue","khaki3","lavender","deeppink2",
          "coral3","beige","brown4","indianred1","lightgreen","orchid")

#' Calculate the mode of vector v
#' @param v A vector.
#' @return The mode.
getmode <- function(v) {

  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Plot the clusters.
#' @param data FIXME
#' @param weights FIXME
#' @param clusters FIXME
#' @param mu FIXME
#' @param main FIXME
#' # FIXME: Consider plotting with ggplot
plot_clusters <- function(data, weights, clusters, mu, main){
  # Plots the clusters

  n <- length(data[,1])
  k <- length(mu[,1])

  # FIXME: Consider plotting with ggplot
  graphics::plot(data, cex = weights / 3,
       pch = 19,
       main = paste(main, " (n = ", n, ")",sep = ""),
       col = c_col[clusters],
       xlim = c(min(data[,1]) - 2, max(data[,1]) + 1),
       ylim = c(min(data[,2]) - 1, max(data[,2]) + 1))
  graphics::points(data[clusters == 99,])
  graphics::points(mu, cex = 2, pch = 4, lwd = 4)
  graphics::legend(min(data[,1]) - 2.5,
                   max(data[,2]) +1,
                   col = c_col[c(1:k, 99)],
                   pch = 19, cex = 0.8,
                   legend = paste("Cluster", c(1:k, 99),
                                  paste("  (", apply(X = t(c(1:k, 99)),
                                                     MARGIN = 2,
                                                     FUN = function(x) {sum(weights[clusters == x])}
                                  ),
                                  ")",
                                  sep=""
                                  )
                   )
  )

}

#' K-means++ algorithm FIXME: insert ref!
#' @param X FIXME
#' @param k FIXME
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

#' Calculate the mode of vector v
#'
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
