#
# Differential expression test for switch like behaviour
# in genes across pseudotime in single-cell RNA-seq data
#
# Convention for naming
# [dist]_[zi?]_[function name]
# e.g. nb_log_likelihood
# or
# e.g. norm_zi_Q_grad etc.
# kieran.campbell@sjc.ox.ac.uk


#' Compute the gradient of the log-likelihood
#' 
#' Computations the gradient of the log-likelihood of the model
#' for params, expression vector x and pseudotime t
#' 
#' @param params A vector with three entries: (1) Parameter \eqn{L = 2 \mu_0}, (2) \eqn{k} 
#' and (3) \code{t_0}
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' 
#' @return The function gradient with respect to \eqn{(L, k, t_0)}.
norm_Q_grad <- function(params, x, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  sig_sq <- params[4]

  ## want to calculate the partial derivatives with respect to the five parameters
  ## L, k, t_0 share the same prefactor vector:
  prefactor <- - 1 / sig_sq * (x - calc_mu(params, t))

  ## L
  dL <- sum(prefactor / (1 + exp(-k*(t - t_0))))

  exp_identity <- 1 / (2 + exp(-k*(t - t_0)) + exp(k*(t - t_0)))

  ## k
  dk <- sum(prefactor * L * (t - t_0) * exp_identity)

  ## t_0
  dt_0 <- sum(-prefactor * L * k * exp_identity)

  ## sigma
  dsig_sq <- -(sum(- 1 / (2 * sig_sq) + 1 / (2 * sig_sq^2) * (x - calc_mu(params, t))^2 ))

  grad <- c(dL, dk, dt_0, dsig_sq)
  if(any(!is.finite(grad))) {
    stop(paste('Infinite gradient at position', which(!is.finite(grad))))
  }
  return(grad)
}


#' Fit sigmoidal expression model
#' 
#' Fit the sigmoidal expression model for some expression vector 
#' x, pseudotime t, and control parameters to be passed to optim
#' This function returns NA if optim fails to converge.
#' 
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' @param control Control arguments to be passed to \code{optim}
norm_fit_alt_model <- function(x, t, control = list(maxit = 100000)) {
  L <- 2 * mean(x) ; t0 <- median(t) ; sig_sq <- var(x)
  k <- coef(lm(x ~ t))[2]
  params <- c(L, k, t0, sig_sq)
  lower <- c(0, -Inf, -Inf, 1e-6)
  opt <- optim(params, fn = norm_alt_obj_func,
               gr = norm_Q_grad,
               method = 'L-BFGS-B',
               lower = lower,
               x = x, t = t, control = control)

  if(opt$convergence == 0) {
    res <- list()
    res$par <- opt$par
    res$par[1] <- res$par[1] / 2
    names(res$par) <- c('mu0','k','t0','sig_sq')
    res$log_lik <- -opt$value
    res$value <- opt$value
    return(res)
  } else {
    warning('Optim failed to converge')
    return( NA )
  }
}

#' Fit null model
#' 
#' Fits the null model (simply by taking the mean and 
#' variance of the expression vector x) and returns a list
#' with the parameters as first entry and negative log likelihood
#' as second.
#' @param x Gene expression vector
norm_fit_null_model <- function(x) {
  mu <- mean(x) ; sig2 <- var(x)
  neg_log_lik <- -log_norm_likelihood(x, mu, sig2)
  par <- c(mu, sig2)
  names(par) <- c("mu", "sig2")
  return(list(par = par, log_lik = -neg_log_lik, value = neg_log_lik))
}

#' Fit sigmoidal and null models
#' 
#' Wrapper to fit both sigmoidal and null models, returning them as a list.
#' 
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' @param ... Additional arguments passed to \code{norm_fit_alt_model}
#' 
#' @return A list of models returned by \code{norm_fit_alt_model} and \code{norm_fit_null_model}
norm_fit_models <- function(x, t, ...) {
  alt_model <- norm_fit_alt_model(x, t, ...)
  if(!is.list(alt_model)) {
    if(!is.na(alt_model)) { # this has to go on a separate line because of how comparisons are done
      stop("alt_model neither list nor NA")
    }
  }
  null_model <- norm_fit_null_model(x)
  return(list(alt_model = alt_model, null_model = null_model))
}

#' Likelihood ratio test for sigmoidal differential expression
#' 
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' @param models List of length two, with the first entry corresponding
#' to the sigmoidal model and the latter to the null model. The model should
#' be of the form of those returned by norm_fit_alt_model and norm_fit_null_model
#' 
#' @return A P-value given by the likelihood ratio test
norm_lrtest <- function(x, t, models) {
  ## first alternative model
  if(length(models$alt_model) < 2) {
    if(is.na(models$alt_model)) return(-1)
  }
  if(length(models$null_model) < 2) {
    if(is.na(models$null_model)) return(-2)
  }
  
  alt_neg_log_lik <- models$alt_model$value
  null_neg_log_lik <- models$null_model$value
  D <- 2 * (null_neg_log_lik - alt_neg_log_lik)
  dof <- 2 # 4 - 2
  return( pchisq(D, dof, lower.tail = FALSE) )
}

#' Sigmoidal differential expression test
#' 
#' Returns the p-value and model for sigmoidal differential expression. 
#' Non-zero-inflated and Gaussian likelihood
#' 
#' @param x Gene expression vector
#' @param t Pseudotime vector
#' @param ... Additional arguments passed to \code{norm_fit_models}
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
#' 
norm_diff_expr_test <- function(x, t, ...) {
  models <- norm_fit_models(x, t, ...)
  pval <- norm_lrtest(x, t, models)
  params <- rep(NA, 4)
  
  if(length(models$alt_model$par == 4)) params <- models$alt_model$par
  
  r <- c(pval, params)
  names(r) <- c('pval', 'mu0', 'k', 't0', 'sig_sq')
  return( r )
}



