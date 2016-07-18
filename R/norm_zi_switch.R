# Expectation maximisation for zero-inflated sigmoidal
# differential expression for single-cell RNA-seq data
# kieran.campbell@sjc.ox.ac.uk

norm_zi_Q <- function(params, x, t, is_zero) {
  par <- params[1:4] ; lambda <- params[5]
  lQ <- alt_obj_func(params, x, t)
  lQ <- lQ + sum(lambda * x[is_zero]^2)
  lQ <- lQ - sum(log(1 - exp(-lambda * x[!is_zero]^2)))
  lQ
}

#' Return the gradient of Q = E[l_c]
norm_zi_Q_grad <- function(params, x, t, is_zero) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  sig_sq <- params[4] ; lambda <- params[5]

  ## want to calculate the partial derivatives with respect to the five parameters
  ## L, k, t_0 share the same prefactor vector:
  prefactor <- - 1 / sig_sq * (x - calc_mu(params, t))

  ## L
  dL <- sum(prefactor / (1 + exp(-k*(t - t_0))))

  ## k
  dk <- sum(prefactor * L * (t - t_0) * exp(-k*(t - t_0)) / (1 + exp(-k*(t - t_0)))^2)

  ## t_0
  dt_0 <- sum(-prefactor * L * k * exp(-k*(t - t_0)) / (1 + exp(-k*(t - t_0)))^2)

  ## sigma
  dsig_sq <- -(sum(- 1 / (2 * sig_sq) + 1 / (2 * sig_sq^2) * (x - calc_mu(params, t))^2 ))

  ## lambda
  x_0 <- x[is_zero] ; x_non_zero <- x[!is_zero]
  dlambda <- -(sum(-x_0^2) +
                 sum(x_non_zero * x_non_zero * exp(-lambda * x_non_zero * x_non_zero) / (1 - exp(-lambda * x_non_zero * x_non_zero))))

  grad <- c(dL, dk, dt_0, dsig_sq, dlambda)
  if(any(!is.finite(grad))) {
    print(c(L, t_0, k, sum(t-t_0), exp(-k*(t - t_0))))
    stop(paste('Infinite gradient at position', which(!is.finite(grad))))
  }
  return(grad)
}


alt_E_step <- function(params, t) {
  lambda <- params[5]
  par <- params[1:4]
  mu <- calc_mu(par, t)
  return( mu / (1 + 2 * par[4] * lambda) )
}

alt_M_step <- function(params, x, t, is_zero, control = list(maxit = 1e6)) {
  ## impose upper bound on t_0
  t_0_max <- 20 * log(10) / params[2]
   opt <- optim(params, fn = norm_zi_Q, gr = norm_zi_Q_grad,
              method = "L-BFGS-B",
              lower = c(0, -Inf, -Inf, 1e-5, 1e-5),
              #upper = c(Inf, Inf, t_0_max, Inf, Inf),
              x = x, t = t, is_zero = is_zero, control = control)
  if(opt$par[3] == t_0_max) warning('t_0 sits at boundary')
  opt
}

#' Runs the EM algorithm for zero inflated model
#' 
#' @param y Gene expression vector
#' @param t Pseudotime vector
#' @param maxit Maximum number of iterations for EM algorithm
#' @param loglik_tol Threshold change in the log likelihood to assume the algorithm has converged.
#' @param verbose Print out EM progress
#' 
#' 
#' @return List with parameters, latent values and log likelihood
alt_EM_ctrl <- function(y, t, maxit = 500, loglik_tol = 1e-7, verbose = FALSE) {
  ## initialisation
  L <- mean(y) ; t_0 <- median(t) ; sig_sq <- var(y)
  k <- coef(lm(y ~ t))[2]
  lambda <- 0.1
  params <- c(L, k, t_0, sig_sq, lambda)

  loglik <- 10^8
  x_latent <- NULL


  for(i in 1:maxit) {
    ## E step
    x <- y
    is_zero <- y == 0
    x[is_zero] <- alt_E_step(params, t)[is_zero]
    x_latent <- x

    ## M step
    opt <- alt_M_step(params, x, t, is_zero)
    
    if(opt$convergence == 0) {
      params <- opt$par
    } else {
      print(opt)
      stop(paste("Something didn't converge:", opt$message))
    }
    delta_loglik <- loglik - opt$value
    if(delta_loglik < 0) stop("negative log likelihood increased - something isn't right")
    if(delta_loglik < loglik_tol) {
      if(verbose) message(paste("EM converged after ", i, " iterations"))
      names(params) <- c('L','k','t_0','sig_sq','lambda')
      return(list(par = params, x = x_latent, loglik = opt$value))
    } else {
      loglik <- opt$value
    }
  }
  warning('EM reached maxit without loglik_tol')
  names(params) <- c('L','k','t_0','sig_sq','lambda')
  return(list(par = params, x = x_latent, loglik = opt$value))
}


# null model EM -----------------------------------------------------------

Null_EM_ctrl <- function(y, loglik_tol = 1e-6, maxit = 100, verbose = FALSE) {
  ## initialisation
  mu <- mean(y) ; sig_sq <- var(y)
  lambda <- 0.1
  params <- c(mu, sig_sq, lambda)

  loglik <- 10^8
  x_latent <- NULL

  for(i in 1:maxit) {
    ## E step
    x <- y
    is_zero <- y == 0
    x[is_zero] <- Null_E_step(params)
    x_latent <- x

    ## M step
    params <- Null_M_step(params, x, is_zero)

    new_loglik <- Null_log_lik(params, x, is_zero)
    delta_loglik <- loglik - new_loglik
    #if(delta_loglik < 0) stop("negative log likelihood increased - something isn't right")
    if(abs(delta_loglik) < loglik_tol) {
      if(verbose) message(paste("EM converged after ", i, " iterations"))
      names(params) <- c('mu','sig_sq','lambda')
      return(list(par = params, x = x_latent, loglik = new_loglik))
    } else {
      loglik <- new_loglik
    }
  }
  warning('EM reached maxit without loglik_tol')
  names(params) <- c('mu','sig_sq','lambda')
  return(list(par = params, x = x_latent, loglik = new_loglik))
}

Null_E_step <- function(params) {
  mu <- params[1] ; sig_sq <- params[2] ; lambda <- params[3]
  return( mu / (1 + 2 * sig_sq * lambda) )
}

lambda_obj_fnc <- function(lambda, x, is_zero) {
    r <- -lambda * sum(x[is_zero])
    r <- r + sum(log(1 - exp(-lambda * x[!is_zero] * x[!is_zero])))
    return( -r )
}


Null_M_step <- function(params, x, is_zero) {
  control = list(maxit = 10000) #, pgtol = 1e-6)
  ## impose upper bound on t_0
  mu <- mean(x) ;
  sig_sq <- 1 / length(x) * sum( (x - mu)^2 )

  ## use old lambda as starting value
  opt <- optim(params[3], fn = lambda_obj_fnc, method = "Brent",
               lower = 0, upper = 1e5,
               x = x, is_zero = is_zero, control = control)
  return(c(mu, sig_sq, opt$par))
}

Null_log_lik <- function(params, x, is_zero) {
  mu <- params[1] ; sig_sq <- params[2] ; lambda <- params[3]
  log_lik <- sum(dnorm(x, mean = mu, sd = sqrt(sig_sq), log = TRUE))
  log_lik <- log_lik + sum(-lambda * x[is_zero])
  log_lik <- log_lik + sum(log(1 - exp(-lambda * x[!is_zero] * x[!is_zero])))
  return(-log_lik)
}

norm_zi_fit_models <- function(x, t, ...) {
  alt_model <- alt_EM_ctrl(x, t, ...)
  null_model <- Null_EM_ctrl(x, ...)
  
  return(list(alt_model = alt_model, null_model = null_model))
}

#' Sigmoidal differential expression test including zero-inflation
#' 
#' Returns the p-value and model for sigmoidal differential expression. 
#' Non-zero-inflated and Gaussian likelihood
#' 
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' 
#' @return A vector of length 5 with entries:
#' \itemize{
#' \item P-value
#' \item MLE estimate for L = 2 mu_0 parameter
#' \item MLE estimate for k
#' \item MLE estimate for t_0
#' \item MLE estimate for sigma^2
#' }
#' 
norm_zi_diff_expr_test <- function(x, t, ...) {
  models <- norm_zi_fit_models(x, t, ...)
  pval <- norm_zi_lrtest(x, t, models)

  params <- models$alt_model$par
  if(length(params) < 5) params <- rep(NA, 5)
  
  r <- c(pval, params)
  names(r) <- c('pval', 'L', 'k', 't0', 'sig_sq', 'lambda')
  return( r )
}

#' Likelihood ratio test for zero-inflated sigmoidal differential expression
#' 
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' @param models List of length two, with the first entry corresponding
#' to the sigmoidal model and the latter to the null model. The model should
#' be of the form of those returned by norm_fit_alt_model and norm_fit_null_model
#' 
#' @return A P-value given by the likelihood ratio test
norm_zi_lrtest <- function(x, t, models) {
  ## first alternative model
  if(length(models$alt_model) < 2) {
    if(is.na(models$alt_model)) return(-1)
  }
  if(length(models$null_model) < 2) {
    if(is.na(models$null_model)) return(-2)
  }
  alt_neg_log_lik <- models$alt_model$loglik
  null_neg_log_lik <- models$null_model$loglik
  D <- 2 * (null_neg_log_lik - alt_neg_log_lik)
  dof <- 2 # 4 - 2
  return( pchisq(D, dof, lower.tail = FALSE) )
}


