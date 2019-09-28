
# Simulated data (n = 200)
# x: x-coordinate
# y: y-coordinate
# w: weights for objects (must be integers)
# id: identification for each object
get_testdata_200 <- function() {
  set.seed(2)
  dat_200 <- data.frame(x = c(stats::rnorm(50, mean = 1.5, sd = 0.8), stats::rnorm(50, mean = 5.5, sd = 0.8),
                              stats::rnorm(50, mean = 3.0, sd = 0.8), stats::rnorm(50, mean = 4.5, sd = 0.8), 8),
                        y = c(stats::rnorm(50, mean = 1.0, sd = 0.8), stats::rnorm(50, mean = 1.5, sd = 1.0),
                              stats::rnorm(50, mean = 4.5, sd= 1.0), stats::rnorm(50, mean = 4.5, sd = 0.8), 6),
                        w =  round(stats::runif(201, min = 0.51, max = 10.49)),
                        id = 1:201)
}

# Simulated data (n = 500)
get_testdata_500 <- function() {
  set.seed(7)
  dat_500 <- data.frame(x = c(stats::rnorm(100, mean = 0.5, sd = 2), stats::rnorm(150, mean = 23.5, sd = 5),
                              stats::rnorm(100, mean = 2.0, sd = 6), stats::rnorm(150, mean = 5.5, sd = 2)),
                        y = c(stats::rnorm(100, mean = 0.5, sd = 2), stats::rnorm(150, mean = 1.0, sd = 7.5),
                              stats::rnorm(100, mean = 20.5, sd= 5.5), stats::rnorm(150, mean = 5.5, sd = 2)),
                        w =  round(stats::runif(500, min = 0.51, max = 10.49)),
                        id = 1:500)
}

# Simulated data (n = 1000)
get_testdata_1000 <- function() {
  set.seed(5)
  dat_1000 <- data.frame(x = c(stats::rnorm(250, mean = 0.5, sd = 2), stats::rnorm(250, mean = 21.5, sd = 3),
                               stats::rnorm(250, mean = 2.0, sd = 6), stats::rnorm(250, mean = 5.5, sd = 2)),
                         y = c(stats::rnorm(250, mean = 0.5, sd = 2), stats::rnorm(250, mean = 1.0, sd = 7.5),
                               stats::rnorm(250, mean = 20.5, sd= 5.5), stats::rnorm(250, mean = 5.5, sd = 2)),
                         w =  round(stats::runif(1000, min = 0.51, max = 10.49)),
                         id = 1:1000)
}

# FIXME: Consider plotting with ggplot
plot_testdata <- function(dat) {
  graphics::plot(dat[, 1:2],
                 cex = dat$w / 3,
                 pch = 19,
                 main = paste("Test data, n =", nrow(dat)),
                 xlim = c(min(dat$x) - 2, max(dat$x) + 1),
                 ylim = c(min(dat$y) - 1 ,max(dat$y) + 1))
}

#' Clustering demo with uniform prior.
#' @export
demo_prob_cluster_uniform <- function() {
  # Number of clusters
  k <- 10

  # Initial mu with k-means
  test_dat <- get_testdata_200()
  init_kmpp <- kmpp(cbind(rep(test_dat$x, test_dat$w),
                          rep(test_dat$y, test_dat$w)),
                    k)
  init_mu <- init_kmpp$centers

  # FIXME: Consider plotting with ggplot
  plot_testdata(test_dat)
  graphics::points(init_mu[, 1], init_mu[, 2], cex = 5, pch = 4, lwd = 5, col = "red")

  # Prior lower und upper limit for cluster size

  # Mean
  pr_mean <- round(sum(test_dat$w) / k)

  # Max width for prior
  pr_width <- 20

  # Lower und upper limit for cluster size
  L <- (pr_mean - pr_width)
  U <- (pr_mean + pr_width)

  # Function call
  lambda1 <- 5
  pt <- proc.time()

  clust <- prob_clust(
    data = test_dat[, 1:2],
    weights = test_dat$w,
    k = k,
    init_mu = init_mu,
    prior_dist = "uniform",
    range = c(L, U),
    lambda = lambda1
  )

  test_dat$cl <- clust[[1]]
  mu <- clust[[2]]
  obj_min <- clust[[3]]
  proc.time() - pt

  plot_clusters(
    test_dat[,1:2],
    test_dat$w,
    test_dat$cl,
    mu,
    title = paste("Probabilistic clustering, uniform distribution on [", L, ", ", U, "], k = ", k, sep = "")
  )
}

#' Clustering demo with Gaussian prior.
#' @export
demo_prob_cluster_normal <- function() {

  # Number of clusters
  k <- 5

  # Initial mu with k-means
  test_dat <- get_testdata_200()
  init_kmpp <- kmpp(cbind(rep(test_dat$x,
                              test_dat$w),
                          rep(test_dat$y,
                              test_dat$w)),
                    k)
  init_mu <- init_kmpp$centers

  # FIXME: Consider plotting with ggplot
  plot_testdata(test_dat)
  graphics::points(init_mu[, 1],
                   init_mu[, 2],
                   cex = 5,
                   pch = 4,
                   lwd = 5,
                   col = "red")

  # Standard deviation
  pr_sd <- 3

  # Function call
  clust <-
    prob_clust(
      data = test_dat[, 1:2],
      weights = test_dat$w,
      k = k,
      init_mu = init_mu,
      prior_dist = "normal",
      sigma = pr_sd
    )
  test_dat$cl <- clust[[1]]
  mu <- clust[[2]]

  # Picture 1
  # FIXME: Consider plotting with ggplot

  plot_clusters(
    test_dat[,1:2],
    test_dat$w,
    test_dat$cl,
    mu,
    title =
      paste("Probabilistic clustering, Gaussian distribution (sd = ", pr_sd, "), k = ", sep = "", k)
  )

}
