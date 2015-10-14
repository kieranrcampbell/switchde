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


#' Log-likelihood of iid measurements from a normal distribution
#' 
#' @param x A vector of measurements
#' @param mu A scalar mean
#' @param sig_sq a scalar variance
#' 
#' @return The log likelihood x given mu, sig_sq
log_norm_likelihood <- function(x, mu, sig_sq) {
  sum(dnorm(x, mu, sqrt(sig_sq), log = TRUE))
}

#' Computations the gradient of the log-likelihood of the model
#' for params, expression vector x and pseudotime t
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


#' Negative log likelihood of the sigmoidal differential expression
#' function given some expression vector x
norm_alt_obj_func <- function(params, x, t) {
  sig_sq <- params[4]
  mu <- calc_mu(params, t)
  return( -log_norm_likelihood(x, mu, sig_sq) )
}

#' Fit the sigmoidal expression model for some expression vector 
#' x, pseudotime t, and control parameters to be passed to optim
#' This function returns NA if optim fails to converge.
norm_fit_alt_model <- function(x, t, control = list(maxit = 100000)) {
  L <- mean(x) ; t0 <- median(t) ; sig_sq <- var(x)
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
    names(res$par) <- c('L','k','t0','sig_sq')
    res$log_lik <- -opt$value
    res$value <- opt$value
    return(res)
  } else {
    warning('Optim failed to converge')
    return( NA )
  }
}

#' Fit's the null model (simply by taking the mean and 
#' variance of the expression vector x) and returns a list
#' with the parameters as first entry and negative log likelihood
#' as second.
norm_fit_null_model <- function(x) {
  mu <- mean(x) ; r <- var(x)
  neg_log_lik <- -log_norm_likelihood(x, mu, r)
  return(list(par = c(mu, r), value = neg_log_lik))
}

#' Fits both sigmoidal and null models, returning them as a list
norm_fit_models <- function(x, t) {
  alt_model <- norm_fit_alt_model(x, t)
  if(!is.list(alt_model) && !is.na(alt_model)) {
    print(alt_model)
    print(class(alt_model))
    stop("alt_model neither list nor NA")
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
#' @export
#' @return A P-value given by the likelihood ratio test
norm_lrtest <- function(x, t, models) {
  ## first alternative model
  if(is.na(models$alt_model)) return(-1)
  if(is.na(models$null_model)) return(-2)
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
#' 
#' @export
#' @return A vector of length 5 with entries:
#' \itemize{
#' \item P-value
#' \item MLE estimate for L = 2 mu_0 parameter
#' \item MLE estimate for k
#' \item MLE estimate for t_0
#' \item MLE estimate for sigma^2
#' }
norm_diff_expr_test <- function(x, t) {
  models <- norm_fit_models(x, t)
  pval <- lrtest(x, t, models)
  params <- rep(NA, 4)
  if(!is.na(models$alt_model)) {
    params <- models$alt_model$par
  }
  r <- c(pval, params)
  names(r) <- c('pval', 'L', 'k', 't0', 'sig_sq')
  return( r )
}
# 
# simulate_de <- function(alt_model, n = 100) {
#   params <- alt_model$par
#   L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
#   t <- runif(n)
#   mu <- L / (1 + exp(-k*(t - t_0)))
#   y <- rnbinom(n, r, mu = mu)
#   return(data.frame(t=t, x=y))
# # }
# 
# simulate_mean <- function(alt_model, n=100) {
#   params <- alt_model$par
#   L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
#   t <- runif(n)
#   mu <- L / (1 + exp(-k*(t - t_0)))
#   return(data.frame(t=t, mu=mu))
# }

norm_plot_model <- function(alt_model, x, t) {
  mu_func <- function(t, params) {
    L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
    L / (1 + exp(-k*(t - t_0)))
  }
  df <- data.frame(y=x, t=t)
  ggplot(df) + geom_point(aes(x=t, y=y), alpha=.8) +
    theme_bw() +
    stat_function(fun = mu_func, args = list(alt_model$par), color='red') +
    xlab('Pseudotime') + ylab('Expression')
}

# testing -----------------------------------------------------------------

test_Q_grad <- function() {
  params <- c(L, k, t_0, sig_sq)
  aof <- function(params, X, t) norm_alt_obj_func(params, x, t)
  grad(aof, x = params, X = X, t = t)
  norm_Q_grad(params, x = x, t = t)
}



