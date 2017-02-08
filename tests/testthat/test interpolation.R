library(SSM)
context("SSM interpolation")

test_that("1d SSM interpolates data", {
  X <- runif(5, -1, 1)
  Y <- runif(5)
  s <- fit.ssm(X, Y)
  predictions <- sapply(X, predict, object = s)
  mapply(expect_equal, Y, predictions)
})

test_that("2d SSM interpolates data", {
  X <- matrix(runif(10, -1, 1), ncol = 2)
  Y <- runif(5)
  s <- fit.ssm(X, Y)
  predictions <- apply(X, 1, predict, object = s)
  mapply(expect_equal, Y, predictions)
})
