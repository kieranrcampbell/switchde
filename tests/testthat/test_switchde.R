
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