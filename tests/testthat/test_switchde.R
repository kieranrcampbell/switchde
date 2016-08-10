
library(switchde)
library(numDeriv)

context("Gradient computations")

test_that("Gradient computation is approximately correct", {
  params <- c(2, 1, 0.5, 0.1)
  n <- 100
  pst <- runif(n)
  x <- rnorm(n, pst^2, params[4])
  x[x < 0] <- 0
  
  analytical_gradient <- grad_likelihood_sigmoidal(params, x, pst)
  numerical_gradient <- grad(likelihood_sigmoidal, x = params, y = x, pst = pst)
  expect_equal(analytical_gradient, numerical_gradient)
})

test_that("EM gradient computation is approximately correct", {
  data(synth_gex)  
  data(ex_pseudotime)
  
  y <- synth_gex[1,]
  pst <- ex_pseudotime
  
  mu0 <- mean(y[y > 0])
  k <- coef(lm(y ~ pst))[2]
  t0 <- median(pst)
  sigma2 <- var(y)
  lambda <- 1
  
  params <- c(mu0, k, t0, sigma2, lambda)

  ## For sigmoidal model
  E <- switchde:::sigmoid_E_step(y, pst, params)
  
  analytical_grad <- switchde:::Q_grad_sigmoid(params, y, pst, E$Ex, E$Ex2)
  numerical_grad <- numDeriv::grad(Q_sigmoid, params, y = y, pst = pst, 
                                   Ex = E$Ex, Ex2 = E$Ex2)
  names(analytical_grad) <- names(numerical_grad) <- NULL
  expect_equal(analytical_grad, numerical_grad)
  
  ## For constant model
  params <- params[c(1,4,5)]
  E <- switchde:::constant_E_step(y, params)
  
  analytical_grad <- switchde:::Q_grad_constant(params, y, E$Ex, E$Ex2)
  numerical_grad <- numDeriv::grad(Q_constant, params, y = y, 
                                   Ex = E$Ex, Ex2 = E$Ex2)
  names(analytical_grad) <- names(numerical_grad) <- NULL
  expect_equal(analytical_grad, numerical_grad)
})



context("Model fitting")

test_that("Model fitting is approximately correct", {
  set.seed(123L)
  n <- 100
  pst <- runif(n)
  k <- 10; t0 <- 0.5; mu0 <- 1
  x <- sigmoid(pst, c(mu0, k, t0))
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

test_that("Output is correct on matrix input for EM algorithm", {
  data(synth_gex)
  data(ex_pseudotime)
  sde <- switchde(synth_gex, ex_pseudotime, zero_inflated = TRUE)
  
  expect_is(sde, "tbl_df")
  expect_is(sde, "data.frame")
  expect_equal(ncol(sde), 7)
  expect_equal(nrow(sde), 12)
})