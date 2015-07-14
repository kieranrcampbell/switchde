#
# Differential expression test for switch like behaviour
# in genes across pseudotime in single-cell RNA-seq data
#
# kieran.campbell@sjc.ox.ac.uk

calc_mu <- function(params, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- L / (1 + exp(-k*(t - t_0)))
  return(mu)
}

#' Log-likelihood of a negative binomial function with count vector k,
#' mean vector mu and constant dispersion parameter r
log_norm_likelihood <- function(x, mu, sig_sq) {
#   N <- length(x)
#   log_lik <- -N / 2 * log(2 * pi * sig_sq)
#   log_lik <- log_lik - sum( 1/(2 * sig_sq) * (x - mu)^2 )
#   return(log_lik)
  sum(dnorm(x, mu, sqrt(sig_sq), log = TRUE))
}

Q_grad <- function(params, x, t) {
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


#' Objective function given differential expression
alt_obj_func <- function(params, x, t) {
  sig_sq <- params[4]
  mu <- calc_mu(params, t)
  return( -log_norm_likelihood(x, mu, sig_sq) )
}

null_obj_func <- function(params, x) {
  mu <- params[1] ; r <- params[2]
  mu_vec <- rep(mu, length(x))
  return( -log_norm_likelihood(x, mu_vec, r) )
}

fit_alt_model <- function(x, t) {
  control <- list(maxit = 10000)
  L <- mean(x) ; t0 <- median(t) ; sig_sq <- var(x)
  k <- coef(lm(x ~ t))[2]
  params <- c(L, k, t0, sig_sq)
  lower <- c(0, -Inf, -Inf, 1e-6)
  opt <- optim(params, fn = alt_obj_func,
               gr = Q_grad,
               method = 'L-BFGS-B',
               lower = lower,
               x = x, t = t, control = control)

  if(opt$convergence == 0) {
    res <- list()
    res$par <- opt$par
    names(res$par) <- c('L','k','t0','sig_sq')
    res$log_lik <- -opt$value
    return(res)
  } else {
    warning('Optim failed to converge')
    return( NA )
  }
}

test_Q_grad <- function() {
#   t <- t0 <- 1
#   x <- 1 ; L <- 2
#   sig_sq <- 1 ; k <- 0 ; lambda <- 1
#   is_zero <- F

  params <- c(L, k, t_0, sig_sq)
  aof <- function(params, X, t) alt_obj_func(params, x, t)
  grad(aof, x = params, X = X, t = t)
  Q_grad(params, x = x, t = t)
}

fit_null_model <- function(x) {
  mu <- mean(x) ; r <- var(x)
  neg_log_lik <- -log_norm_likelihood(x, mu, r)
  return(list(par = c(mu, r), value = neg_log_lik))
}

fit_models <- function(x, t) {
  alt_model <- fit_alt_model(x, t)
  null_model <- fit_null_model(x)
  return(list(alt_model = alt_model, null_model = null_model))
}

lrtest <- function(x, t, models) {
  ## first alternative model
#   if(is.na(models$alt_model)) return(1)
#   if(is.na(models$null_model)) return(0)
  alt_neg_log_lik <- models$alt_model$value
  null_neg_log_lik <- models$null_model$value
  D <- 2 * (null_neg_log_lik - alt_neg_log_lik)
  dof <- 2 # 4 - 2
  return( pchisq(D, dof, lower.tail = FALSE) )
}

diff_expr_test <- function(x, t) {
  models <- fit_models(x, t)
  pval <- lrtest(x, t, models)
  params <- models$alt_model$par
  r <- c(pval, params)
  names(r) <- c('pval', 'L', 'k', 't0', 'sig_sq')
  return( r )
}

simulate_de <- function(alt_model, n = 100) {
  params <- alt_model$par
  L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
  t <- runif(n)
  mu <- L / (1 + exp(-k*(t - t_0)))
  y <- rnbinom(n, r, mu = mu)
  return(data.frame(t=t, x=y))
}

simulate_mean <- function(alt_model, n=100) {
  params <- alt_model$par
  L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
  t <- runif(n)
  mu <- L / (1 + exp(-k*(t - t_0)))
  return(data.frame(t=t, mu=mu))
}

plot_model <- function(alt_model, x, t) {
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


