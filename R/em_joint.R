# Expectation maximisation for zero-inflated sigmoidal
# differential expression for single-cell RNA-seq data
# kieran.campbell@sjc.ox.ac.uk


# differential expression EM ----------------------------------------------

#' Taken from package 'kyotil'
logDiffExp <- function(logx1, logx2, c=1) {
  if (logx1 < logx2) {
    cat("\nlWarning [logDiffExp]: first argument smaller than second, return NaN.\n\n")
    return(NaN)
  }
  #c = max(logx1, logx2)
  log(exp(logx1 - c) - exp(logx2 - c)) + c
}

to_matrix <- function(params) {
  if(length(params) %% 4 != 1) stop('params must be of length 4n+1')
  np <- length(params) - 1
  return(list(matrix(params[1:np], ncol = 4, byrow = TRUE), params[length(params)]))
}

to_vector <- function(params) {
  return(c(t(params[[1]]), params[[2]]))
}

calc_mu <- function(params, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- L / (1 + exp(-k*(t - t_0)))
}

log_norm_likelihood <- function(x, mu, sig_sq) {
#   N <- length(x)
#   log_lik <- -N / 2 * log(2 * pi * sig_sq)
#   log_lik <- log_lik - sum( 1/(2 * sig_sq) * (x - mu)^2 )
#   return(log_lik)
  sum(dnorm(x, mu, sqrt(sig_sq), log = TRUE))
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
  log_sum_expr <- log(1 - exp(-lambda * x[!is_zero]^2))
  is_undef <- !is.finite(log_sum_expr)
  if(any(is_undef))  log_sum_expr[is_undef] <- log(lambda * x[!is_zero][is_undef]^2)
  lQ <- lQ - sum(log_sum_expr)
  if(!is.finite(lQ)) stop('About to return non-finite E[l_c]')
  lQ
}

Q_joint <- function(params, X, t, is_zero) {
  p_mat <- to_matrix(params)
  P <- do.call(cbind, p_mat) # add lambda as final parameter
  Qs <- sapply(1:nrow(P), function(i) {
    Q(P[i,], X[i,], t, is_zero[i,])
    #print(i)
  })
  Qj <- sum(Qs)
  if(!is.finite(Qj)) stop('Trying to return non-finite Q')
  return( Qj )
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

  exp_identity <- 1 / (2 + exp(-k*(t - t_0)) + exp(k*(t - t_0)))

  ## k
  dk <- sum(prefactor * L * (t - t_0) * exp_identity)

  ## t_0
  dt_0 <- sum(-prefactor * L * k * exp_identity)

  ## sigma
  dsig_sq <- -(sum(- 1 / (2 * sig_sq) + 1 / (2 * sig_sq^2) * (x - calc_mu(params, t))^2 ))

  grad <- c(dL, dk, dt_0, dsig_sq)
  if(any(!is.finite(grad))) {
    print(c(L, t_0, k, sum(t-t_0), exp(-k*(t - t_0))))
    stop(paste('Infinite gradient at position', which(!is.finite(grad))))
  }
  return(grad)
}

complete_lambda_grad <- function(X, lambda, is_zero) {
  x_0 <- as.vector(X[is_zero])
  x_n0 <- as.vector(X[!is_zero])
  g <- sum(-x_0 * x_0)
  second_term <- x_n0^2 / (exp(lambda * x_n0^2) - 1)
  if(any(!is.finite(second_term))) {
    undef <- which(!is.finite(second_term))
    second_term[undef] <- 1 / lambda - x_n0[undef]^2 / 2
  }
  g <- g + sum(second_term)
  if(!is.finite(g)) stop('Trying to return non-finite lambda derivative')
  return(-g)
}

Q_grad_joint <- function(params, X, t, is_zero) {
  p_mat <- to_matrix(params)
  P <- do.call(cbind, p_mat) # add lambda as final parameter
  Qs <- c(sapply(1:(nrow(P)), function(i) {
    Q_grad(P[i,], X[i,], t, is_zero[i,])
  }))
  grad <- c(Qs, complete_lambda_grad(X, p_mat[[2]], is_zero))
  grad
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

test_Q_joint_grad <- function() {
#   t <- t0 <- 1
#   x <- 1 ; L <- 2
#   sig_sq <- 1 ; k <- 0 ; lambda <- 1
#   params <- c(L, k, t0, sig_sq, lambda)
  L <- rowMeans(Y) ; t_0 <- median(t) ; sig_sq <- apply(Y, 1, var)
  k <- apply(Y, 1, function(y) coef(lm(y ~ t))[2])
  lambda <- 4

  params <- cbind(L, k, rep(t_0, nrow(Y)), sig_sq)
  params <- to_vector(list(params, lambda))

  #is_zero <- F
  Q_joint(params, X=X, t = t, is_zero = is_zero)
  grad(Q_joint, x = params, X = X, t = t, is_zero = is_zero)
  Q_grad_joint(params, X = X, t = t, is_zero = is_zero)
}


E_step <- function(params, t) {
  p_mat <- to_matrix(params)
  lambda <- p_mat[[2]]
  Mu <- apply(p_mat[[1]], 1, calc_mu, t) # pst-by-gene
  sigmas <- p_mat[[1]][,4]
  apply(Mu, 1, function(mu) mu / (1 + 2 * lambda * sigmas)) # gene-by-pst
}

M_step <- function(params, X, t, is_zero) {
  control = list(maxit = 20000)#, factr = 1e5, pgtol = 1e-6)
  ## impose upper bound on t_0
  t_0_max <- 20 * log(10) / params[2]
  lower <-  c(rep(c(0, -Inf, -Inf, 1e-5), nrow(X)), 1e-5)
  ## upper <- c(rep(c(Inf, Inf, t_0_max, Inf), nrow(X)), Inf) # not needed after numerical correction
  opt <- optim(params, fn = Q_joint, gr = Q_grad_joint,
              method = "L-BFGS-B",
              lower = lower,
              #upper = upper,
              X = X, t = t, is_zero = is_zero, control = control)
  ## if(opt$par[3] == t_0_max) warning('t_0 sits at boundary')
  opt
}

EM_ctrl <- function(Y, t, loglik_tol = 1e-6, maxit = 100) {
  ## initialisation
  L <- rowMeans(Y) ; t_0 <- median(t) ; sig_sq <- apply(Y, 1, var)
  k <- apply(Y, 1, function(y) coef(lm(y ~ t))[2])
  lambda <- 1

  params <- cbind(L, k, rep(t_0, nrow(Y)), sig_sq)
  params <- to_vector(list(params, lambda))

  loglik <- 10^20
  x_latent <- NULL

  for(i in 1:maxit) {
    ## print(paste('EM step', i))

    ## E step
    X <- Y
    is_zero <- Y == 0
    X[is_zero] <- E_step(params, t)[is_zero]
    #x_latent <- x

    ## M step
    opt <- M_step(params, X, t, is_zero)
    if(opt$convergence == 0) {
      params <- opt$par
    } else {
      print(opt)
      numeric_grad <- grad(Q_joint, x = params, X = X, t = t, is_zero = is_zero)
      analytic_grad <- Q_grad_joint(params, X, t, is_zero)
      ## print(qplot(numeric_grad, analytic_grad))
      ## plot(numeric_grad - analytic_grad)
      stop(paste("Something didn't converge:", opt$message))
    }
    delta_loglik <- loglik - opt$value
    if(delta_loglik < 0) stop("negative log likelihood increased - something isn't right")
    if(delta_loglik < loglik_tol) {
      #names(params) <- c('L','k','t_0','sig_sq','lambda')
      message(paste('Completed expectation-maximisation in', i, 'iterations'))
      return(list(par = to_matrix(params), loglik = opt$value))
    } else {
      loglik <- opt$value
    }
  }
  warning('EM reached maxit without loglik_tol')
  names(params) <- c('L','k','t_0','sig_sq','lambda')
  return(list(par = to_matrix(params), loglik = opt$value))
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

plot_joint_models <- function(params, Y, t) {
  library(gridExtra)
  library(ggplot2)
  n_genes <- nrow(Y)
  plots <- lapply(1:n_genes, function(i) {
    y <- Y[i,]
    par <- c(params[[1]][i,], params[[2]])
    plt <- plot_model(par, y, t) + ggtitle(rownames(Y)[i])
    plt
  })
  do.call(grid.arrange, plots)
}


