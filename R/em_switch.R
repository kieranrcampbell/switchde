# Expectation maximisation for zero-inflated sigmoidal
# differential expression for single-cell RNA-seq data
# kieran.campbell@sjc.ox.ac.uk


# differential expression EM ----------------------------------------------

calc_mu <- function(params, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- L / (1 + exp(-k*(t - t_0)))
  return( mu )
}

log_norm_likelihood <- function(x, mu, sig_sq) {
  N <- length(x)
  log_lik <- -N / 2 * log(2 * sig_sq)
  log_lik <- log_lik - sum( 1/(2 * sig_sq) * (x - mu)^2 )
  log_lik <- log_lik - N / 2 * log(2 * pi)
  return(log_lik)
}

#' Objective function given differential expression
alt_obj_func <- function(params, x, t) {
  mu <- calc_mu(params, t)
  sig_sq <- params[4]
  return( -log_norm_likelihood(x, mu, sig_sq) )
}

Q <- function(params, x, t, is_zero) {
  par <- params[1:4] ; lambda <- params[5]
  lQ <- alt_obj_func(params, x, t)
  lQ <- lQ + sum(lambda * x[is_zero]^2)
  lQ <- lQ - sum(log(1 - exp(-lambda * x[!is_zero]^2)))
  lQ
}

#' Return the gradient of Q = E[l_c]
Q_grad <- function(params, x, t, is_zero) {
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

test_Q_grad <- function() {
  t <- t0 <- 1
  x <- 1 ; L <- 2
  sig_sq <- 1 ; k <- 0 ; lambda <- 1
  params <- c(L, k, t0, sig_sq, lambda)
  is_zero <- F
  fdHess(params, Q, x = x, t = t, is_zero = is_zero)$gradient
  Q_grad(params, x = x, t = t, is_zero = is_zero)
}

E_step <- function(params, t) {
  lambda <- params[5]
  par <- params[1:4]
  mu <- calc_mu(par, t)
  return( mu / (1 + 2 * par[4] * lambda) )
}

M_step <- function(params, x, t, is_zero) {
  control = list(maxit = 10000) #, pgtol = 1e-6)
  ## impose upper bound on t_0
  t_0_max <- 20 * log(10) / params[2]
   opt <- optim(params, fn = Q, gr = Q_grad,
              method = "L-BFGS-B",
              lower = c(0, -Inf, -Inf, 1e-5, 1e-5),
              #upper = c(Inf, Inf, t_0_max, Inf, Inf),
              x = x, t = t, is_zero = is_zero, control = control)
  if(opt$par[3] == t_0_max) warning('t_0 sits at boundary')
  opt
}

EM_ctrl <- function(y, t, loglik_tol = 1e-6, maxit = 100) {
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
    x[is_zero] <- E_step(params, t)[is_zero]
    x_latent <- x

    ## M step
    opt <- M_step(params, x, t, is_zero)
    if(opt$convergence == 0) {
      params <- opt$par
    } else {
      print(opt)
      stop(paste("Something didn't converge:", opt$message))
    }
    delta_loglik <- loglik - opt$value
    if(delta_loglik < 0) stop("negative log likelihood increased - something isn't right")
    if(delta_loglik < loglik_tol) {
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

Null_EM_ctrl <- function(y,loglik_tol = 1e-6, maxit = 100) {
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
    #print(paste("New negative log likelihoo:", new_loglik))
    #print(params)
    delta_loglik <- loglik - new_loglik
    #if(delta_loglik < 0) stop("negative log likelihood increased - something isn't right")
    if(abs(delta_loglik) < loglik_tol) {
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

plot_model <- function(params, x, t) {
  mu_func <- function(t, params) {
    L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
    L / (1 + exp(-k*(t - t_0)))
  }
  df <- data.frame(y=x, t=t)
  ggplot(df) + geom_point(aes(x=t, y=y), alpha=.8) +
    theme_bw() +
    stat_function(fun = mu_func, args = list(params), color='red') +
    xlab('Pseudotime') + ylab('Expression')
}

plot_null_model <- function(params, x, t) {
  mu <- params[1]
  df <- data.frame(y=x, t=t, mu = mu)
  ggplot(df) + geom_point(aes(x=t, y=y), alpha=.8) +
    theme_bw() +
    geom_line(aes(x=t, y=mu), color='red') +
    xlab('Pseudotime') + ylab('Expression')
}

# library(monocle)
# data(HSMM)
# cds <- HSMM[, HSMM$State %in% 1:2]
#
# mrf_genes <- c('CDK1', 'MEF2C', 'MYH3', 'MYOG', 'ID1')
# mrf_long_genes <- row.names(fData(HSMM))[match(mrf_genes, fData(HSMM)$gene_short_name)]
# mrf_indices <- match(mrf_long_genes, featureNames(cds))
# mrf_expr <- exprs(cds)[mrf_indices,]
#
pst <- cds$Pseudotime
y <- log2(mrf_expr[1,] + 1)
qplot(pst, y)
EM <- EM_ctrl(y, pst, loglik_tol = 1e-7)
plot_model(EM$par, y[s], pst[s])

Null_EM <- Null_EM_ctrl(y)
plot_null_model(Null_EM$par, y, pst)

D <- 2*(Null_EM$loglik - EM$loglik)
pchisq(D, 2, lower.tail = F)

# qplot(pst, EM$x)
#
# i <- 1
# plots <- apply(mrf_expr, 1, function(e) {
#   y <- log2(e+1)
#   EM <- EM_ctrl(y, pst)
#   print(EM$par[2])
#   plt <- plot_model(EM$par, EM$x, pst) + ggtitle(mrf_genes[i]) #+ scale_y_log10()
#   i <<- i + 1
#   plt
# })
#
# library(gridExtra)
# do.call(grid.arrange, plots)
#
# # bootstrap k -------------------------------------------------------------
df <- data.frame(y = log2(mrf_expr[1,] + 1), pst <- cds$Pseudotime)
stat_wrapper <- function(data, indices)  {
  df <- data[indices,]
  EM_ctrl(df$y, df$pst)$par[2]
}
bt <- jackknife(df, stat_wrapper, 500)

