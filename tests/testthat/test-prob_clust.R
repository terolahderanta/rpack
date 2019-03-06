context("test-prob_clust")

test_that("arguments are properly checked", {
  expect_error(prob_clust(NULL, NULL, NULL))
  expect_error(prob_clust("foo", NULL, NULL))
  expect_error(prob_clust(tibble(), "foo", "foo"))
  expect_error(prob_clust(tibble(), c(), "foo"))
  expect_error(prob_clust(tibble(x1 = c(0)), c(1), "foo"))
  expect_error(prob_clust(tibble(x1 = c(0)), c(1), 3))
  expect_error(prob_clust(tibble(x1 = c(0, 1)), c(1, 2), 3))
  expect_error(prob_clust(tibble(x1 = c(0, 1)), c(1, 2), 2))

  # These should work?
  expect_silent(prob_clust(data = tibble::tibble(x1 = c(0, 1), x2 = c(1, 2)),
                           weights = c(1, 2),
                           k = 2))
  expect_silent(prob_clust(data = tibble::tibble(x1 = c(0, 1, 2), x2 = c(1, 2, 1)),
                           weights = c(1, 2, 1),
                           k = 2))
})

