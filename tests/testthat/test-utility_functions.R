context("test-utility_functions")

test_that("getmode works", {
  expect_equal(getmode(1), 1)
  expect_equal(getmode(c(1,1)), 1)
  expect_equal(getmode(c(1,1,1)), 1)
  expect_equal(getmode(c(2,1,1)), 1)
  expect_equal(getmode(c(1,2,1)), 1)
  expect_equal(getmode(c(1,1,2)), 1)
  expect_equal(getmode(c(2,1,3,5,2,4,2,3)), 2)
})

test_that("kmpp cluster initialization works", {
  set.seed(2)

  dat <- tibble::tibble(
    x1 = c(4, 5, 1, 2, 8),
    x2 = c(1, 3, 2, 5, 3)
  )
  w <-  c(1, 1, 3, 2, 1)
  k <-  3

  initial_centers <- kmpp(dat, k)$initial_centers

  expect_equal(initial_centers,
               matrix(
                 c(4, 1,
                   1, 2,
                   2, 5),
                 nrow = 3, ncol = 2,
                 byrow = TRUE,
                 dimnames = list(NULL,c("x1","x2"))
               )
  )
})
