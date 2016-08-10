## Zero inflated sigmoidal differential expression
## kieran.campbell@sjc.ox.ac.uk


# Conventions:
# 
# We consider 
# - data vectors y and pst (observed gene counts and pseudotimes)
# - latent variable x
# - parameters [mu0, k, t0, sigma2, lambda] - *always in this order *

#' @importFrom stats coef median var lm optim
EM_sigmoid <- function(y, pst, iter = 100, log_lik_tol = 1e-3,
                       verbose = FALSE) {
  ## Initialise
  # pars <- initialise()
  mu0 <- mean(y[y > 0])
  k <- coef(lm(y ~ pst))[2]
  t0 <- median(pst)
  sigma2 <- var(y)
  lambda <- 1
  
  params <- c(mu0, k, t0, sigma2, lambda)
  
  param_bounds <- c(0.1, -Inf, -Inf, 0.1, 0.01)
  
  Q_val <- -Inf
  
  ## EM algorithm
  for(it in seq_len(iter)) {
    # E step
    E <- sigmoid_E_step(y, pst, params)
    
    
    # M step
    opt <- optim(params, 
                 fn = Q_sigmoid, gr = Q_grad_sigmoid, 
                 y = y, pst = pst, Ex = E$Ex, Ex2 = E$Ex2,
                 method = "L-BFGS-B", lower = param_bounds, 
                 control = list(fnscale = -1))
    params <- opt$par
    if(abs(Q_val - opt$value) < log_lik_tol) {
      # Converged
      if(verbose) message(paste("Expectation-maximisation converged in", it, "iterations"))
      return(list(params = params, x = E$Ex, log_lik = opt$value))
    }
    Q_val <- opt$value
  }
  warning("EM algorithm failed to converge. Consider increasing maximum iterations.")
  warning("Returning most recent parameter estimates anyway")
  return(list(params = opt$params, x = E$Ex, log_lik = opt$value))
}

sigmoid_E_step <- function(y, pst, params) {
  mu <- sigmoid(pst, params)
  is_zero <- y == 0
  
  Ex <- y
  Ex2 <- y^2
  
  x_expectation <- function(mu, sigma2, lambda) {
    return( mu / (2 * sigma2 * lambda + 1))
  }
  x2_expectation <- function(mu, sigma2, lambda) {
    return(rep(sigma2 / (2 * sigma2 * lambda + 1), length(mu)))
  }
  
  Ex[is_zero] <- sapply(mu[is_zero], x_expectation, params[4], params[5])
  Ex2[is_zero] <- sapply(mu[is_zero], x2_expectation, params[4], params[5])
  
  return(list(Ex = Ex, Ex2 = Ex2))
}



# Taken from package 'kyotil'
logDiffExp <- function(logx1, logx2, c=1) {
  log(exp(logx1 - c) - exp(logx2 - c)) + c
}

Q_sigmoid <- function(params, y, pst, Ex, Ex2) {
  # Setup
  is_zero <- y == 0
  sigma2 <- params[4]; lambda <- params[5]
  mu <- sigmoid(pst, params)
  C <- length(mu)  
  
  # Begin Q calc
  Q <- - C / 2 * log(2 * pi * sigma2)  
  qtmp <- - (1 / (2 * sigma2) + lambda) * Ex2[is_zero] +
    mu[is_zero] * Ex[is_zero] / sigma2 -
    mu[is_zero]^2 / (2 * sigma2)
  Q <- Q + sum(qtmp)

  qtmp <- -(y[!is_zero] - mu[!is_zero])^2 / (2 * sigma2) +
    logDiffExp(0, -lambda * y[!is_zero]^2)
  Q <- Q + sum(qtmp)
  
  if(!is.finite(Q)) print(params)

  return(Q)
}

Q_grad_sigmoid <- function(params, y, pst, Ex, Ex2) {
  mu0 <- params[1]; k <- params[2]; t0 <- params[3]
  sigma2 <- params[4]; lambda <- params[5]
  C <- length(y)
  is_zero <- y == 0
  
  mu <- sigmoid(pst, params)
  
  dQdmu <- 1 / sigma2 * (Ex - mu) # Ex = y if y > 0
  
  dmudmu0 <- mu / mu0
  dmudk <- mu * (pst - t0) / (1 + exp(k * (pst - t0)))
  dmudt0 <- - mu * k / (1 + exp(k * (pst - t0)))
  
  dQdmu0 <- sum(dQdmu * dmudmu0)
  dQdk <- sum(dQdmu * dmudk)
  dQdt0 <- sum(dQdmu * dmudt0)

  dQdsigma2 <- -C / (2 * sigma2) + 
    1 / (2 * sigma2^2) * sum(Ex2 - 2 * mu * Ex + mu^2)

  dQdlambda <- -sum(Ex2[is_zero]) +
    sum(y[!is_zero]^2 / (exp(lambda * y[!is_zero]^2) - 1))

  return(c(dQdmu0, dQdk, dQdt0, dQdsigma2, dQdlambda))  
}

