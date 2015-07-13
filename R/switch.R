# 
# Differential expression test for switch like behaviour
# in genes across pseudotime in single-cell RNA-seq data
# 
# kieran.campbell@sjc.ox.ac.uk

#' Log-likelihood of a negative binomial function with count vector k,
#' mean vector mu and constant dispersion parameter r
log_neg_bin_likelihood <- function(k, mu, r) {
  N <- length(k)
  k_sum <- sum(k)
  log_lik <- N * (r * log(r) - lgamma(r))
  log_lik <- log_lik + sum(lgamma(r + k))
  log_lik <- log_lik - sum(r * log(r + mu))
  log_lik <- log_lik + sum(k * log(mu))
  log_lik <- log_lik - sum(k * log(r + mu))
  log_lik
}

#' Objective function given differential expression
alt_obj_func <- function(params, x, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
  mu <- L / (1 + exp(-k*(t - t_0)))
  return( -log_neg_bin_likelihood(x, mu, r) )
}

null_obj_func <- function(params, x) {
  mu <- params[1] ; r <- params[2]
  mu_vec <- rep(mu, length(x))
  return( -log_neg_bin_likelihood(x, mu_vec, r) )
}

fit_alt_model <- function(x, t) {
  control <- list(maxit = 10000)
  L <- mean(x) ; t_0 <- median(t) ; r <- 1
  k <- coef(lm(x ~ t))[2]
  opt <- optim(c(L, k, t_0, r), fn = alt_obj_func, 
               x = x, t = t, control = control)
  if(opt$convergence == 0) {
    return(opt)
  } else {
    warning('Optim failed to converge')
    return( NA )
  }
}

fit_null_model <- function(x) {
  mu <- mean(x) ; r <- 1
  opt <- optim(c(mu, r), fn = null_obj_func, x = x)
  if(opt$convergence == 0) {
    return(opt)
  } else {
    warning('Optim failed to converge')
    return( NA )
  }
}

fit_models <- function(x, t) {
  alt_model <- fit_alt_model(x, t)
  null_model <- fit_null_model(x)
  return(list(alt_model = alt_model, null_model = null_model))
}

lrtest <- function(x, t, models) {
  ## first alternative model
  if(is.na(models$alt_model)) return(1)
  if(is.na(models$null_model)) return(0)
  alt_neg_log_lik <- models$alt_model$value
  null_neg_log_lik <- models$null_model$value
  D <- 2 * (null_neg_log_lik - alt_neg_log_lik)
  dof <- 2 # 4 - 2
  return( pchisq(D, dof, lower.tail = FALSE))
}

diff_expr_test <- function(x, t) {
  models <- fit_models(x, t)
  lrtest(x, t, models)
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


