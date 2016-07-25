
# Maximum likelihood estimates for non zero-inflated model
# As before params = (mu0, k, t0, sigma2)

#' @importFrom stats dnorm
likelihood_sigmoidal <- function(params, y, pst) {
  mu <- sigmoid(pst, params[1:3])
  return(sum(dnorm(y, mu, sqrt(params[4]), log = TRUE)))
}

grad_likelihood_sigmoidal <- function(params, y, pst) {
  C <- length(y)
  mu <- sigmoid(pst, params[1:3])  
  mu0 <- params[1]; k <- params[2]; 
  t0 <- params[3]; sigma2 <- params[4]
  
  dmudmu0 <- mu / mu0
  dmudk <- mu * (pst - t0) / (1 + exp(k * (pst - t0)))
  dmudt0 <- - mu * k / (1 + exp(k * (pst - t0)))
  
  pre_vec <- 1 / params[4] * (y - mu)
  
  dfdmu0 <- sum(pre_vec * dmudmu0)
  dfdk <- sum(pre_vec * dmudk)
  dfdt0 <- sum(pre_vec * dmudt0)
  
  dfdsigma2 <- - C / (2 * sigma2) + 1 / (2 * sigma2^2) * sum((y - mu)^2)

  return(c(dfdmu0, dfdk, dfdt0, dfdsigma2))
}

likelihood_constant <- function(params, y) {
    mu <- params[1]; sigma2 <- params[2]
    return(sum(dnorm(y, mu, sqrt(sigma2), log = TRUE)))
}

fit_sigmoidal_model <- function(y, pst) {
    mu0 <- mean(y[y > 0])
    k <- coef(lm(y ~ pst))[2]
    t0 <- median(pst)
    sigma2 <- var(y)
    
    params <- c(mu0, k, t0, sigma2)
    param_bounds <- c(0.1, -Inf, -Inf, 0.1)
    
    opt <- optim(params, 
                 fn = likelihood_sigmoidal, gr = grad_likelihood_sigmoidal, 
                 y = y, pst = pst, 
                 method = "L-BFGS-B", lower = param_bounds, 
                 control = list(fnscale = -1))
    
    params <- opt$par
    names(params) <- c("mu0", "k", "t0", "sigma2")
    return(list(params = params, log_lik = opt$value))
}

fit_constant_model <- function(y) {
  mu <- mean(y)
  sigma2 <- var(y)
  params <- c(mu, sigma2)
  names(params) <- c("mu", "sigma2")
  log_lik <- likelihood_constant(params, y)
  return(list(params = params, log_lik = log_lik))
}

