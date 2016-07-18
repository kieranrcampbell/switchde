
library(switchde)
library(numDeriv)

context("Gradient computations")

test_that("Gradient computation is approximately correct", {
  params <- c(2, 1, 0.5, 1)
  n <- 100
  x <- rnorm(n, 0, 1)
  x[x < 0] <- 0
  pst <- runif(n)
  
  analytical_gradient <- switchde:::norm_Q_grad(params, x = x, t = pst)
  aof <- function(params, xx, t) switchde:::norm_alt_obj_func(params, xx, t)
  numerical_gradient <- grad(aof, x = params, xx = x, t = pst)
  expect_equal(analytical_gradient, numerical_gradient)
})

context("Model fitting")

test_that("Model fitting is approximately correct", {
  set.seed(123L)
  n <- 100
  pst <- runif(n)
  k <- 10; t0 <- 0.5; mu0 <- 1
  x <- switchde:::calc_mu(c(2 * mu0, k, t0), pst)
  x <- rnorm(n, x, 0.1)
  x[x < 0] <- 0
  sde <- switchde(x, pst)
  
  ## Check sde is correct class
  expect_is(sde, "tbl_df")
  expect_is(sde, "data.frame")

  expect_equal(ncol(sde), 6)
  
  ## Check parameter estimates are roughly correct &
  ## p-val is less that 0.05
  expect_lt(sde$pval, 0.05)
  
  ## check parameter estimates aren't wildly off
  expect_lt(abs(sde$mu0 - mu0), 0.1)
  expect_lt(abs(sde$k - k), 1) # more tolerance for k due to geometry
  expect_lt(abs(sde$t0 - t0), 0.1)
})

test_that("Output is correct on matrix input", {
  data(example_gex)
  data(example_pseudotime)
  sde <- switchde(example_gex, example_pseudotime)

  expect_is(sde, "tbl_df")
  expect_is(sde, "data.frame")
  expect_equal(ncol(sde), 6)
  expect_equal(nrow(sde), 4)
})