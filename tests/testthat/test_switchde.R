
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
})

test_that("Output is correct on matrix input", {
  data(synth_gex)
  data(ex_pseudotime)
  sde <- switchde(synth_gex, ex_pseudotime)

  expect_is(sde, "tbl_df")
  expect_is(sde, "data.frame")
  expect_equal(ncol(sde), 6)
  expect_equal(nrow(sde), 12)
})