library(SSM)
context("Computing variances using Q")

test_that("the Q-computed variances match the sum of term variances", {
  X <- matrix(runif(10, -1, 1), ncol = 2)
  Y <- runif(5)
  s <- fit.ssm(X, Y, SA = TRUE)
  Q <- constructQ(s,1)
  expect_equal(as.numeric(t(Y) %*% Q %*% Y), s@main_sobol[1] * s@total_variance)
})
