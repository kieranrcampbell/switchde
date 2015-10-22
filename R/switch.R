#
# Differential expression test for switch like behaviour
# in genes across pseudotime in single-cell RNA-seq data
#
# kieran.campbell@sjc.ox.ac.uk

#' Log-likelihood of a negative binomial function with count vector k,
#' mean vector mu and constant dispersion parameter r
log_neg_bin_likelihood <- function(k, mu, r) {
  return( sum(dnbinom(k, size = r, mu = mu, log = TRUE)) )
}

log_norm_likelihood <- function(x, mu, sig_sq) {
  sum(dnorm(x, mu, sqrt(sig_sq), log = TRUE))
}

calc_mu <- function(params, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- L / (1 + exp(-k*(t - t_0)))
  return(mu)
}

#' Objective function given differential expression
nb_alt_obj_func <- function(params, x, t) {
  r <- params[4]
  mu <- calc_mu(params, t)
  return( -log_neg_bin_likelihood(x, mu, r) )
}

nb_null_obj_func <- function(params, x) {
  mu <- params[1] ; r <- params[2]
  mu_vec <- rep(mu, length(x))
  return( -log_neg_bin_likelihood(x, mu_vec, r) )
}

#' Objective function given differential expression
norm_alt_obj_func <- function(params, x, t) {
  sig_sq <- params[4]
  mu <- calc_mu(params, t)
  return( -log_norm_likelihood(x, mu, sig_sq) )
}

norm_null_obj_func <- function(params, x) {
  mu <- params[1] ; r <- params[2]
  mu_vec <- rep(mu, length(x))
  return( -log_norm_likelihood(x, mu_vec, r) )
}

fit_alt_model <- function(x, t, obj_func, dist = c('nb','normal')) {
  if(!is.integer(x) && (dist == 'nb')) {
    warning('Fitting NB distribution with non-integer counts. Rounding x')
    x <- round(x)
  }
  control <- list(maxit = 10000)
  L <- mean(x) ; t_0 <- median(t) ; r <- 1
  k <- coef(lm(x ~ t))[2]
  params <- c(L, k, t_0, r)
  opt <- tryCatch({
    optim(params, fn = obj_func,
               x = x, t = t, control = control)
  }, error = function(e) {
    print(paste('Optimisation failed:', e))
    return( NULL )
  })
  if(is.null(opt)) return(NULL)
  if(opt$convergence == 0) {
    res <- list()
    res$par <- opt$par
    names(res$par) <- c('L','k','t0',
                        ifelse(dist == 'nb', 'r', 'sig_sq'))
    res$log_lik <- -opt$value
    res$df <- 4
    return(res)
  } else {
    warning('Optim failed to converge')
    return( NA )
  }
}

fit_null_model <- function(x, obj_func, dist = c('nb','normal')) {
  if(!is.integer(x) & (dist == 'nb')) {
    warning('Fitting NB distribution with non-integer counts. Rounding x')
    x <- round(x)
  }
  mu <- mean(x) ; r <- 1
  opt <- tryCatch({
    optim(c(mu, r), fn = obj_func, x = x)
  }, error = function(e) {
    print(paste('Optimisation failed:', e))
    return( NULL )
  })
  if(is.null(opt)) return(NULL)
  if(opt$convergence == 0) {
    res <- list()
    res$par <- opt$par
    names(res$par) <- c('mu',ifelse(dist == 'nb', 'r', 'sig_sq'))
    res$log_lik <- -opt$value
    res$df <- 2
    return(res)
  } else {
    warning('Optim failed to converge')
    return( NULL )
  }
}

# #' Fit pseudotime model
# #'
# #' @export
# fitModel <- function(x, t = NULL, dist = c('nb', 'normal'), sigmoidal = TRUE) {
#   dist <- match.arg(dist)
#   if(sigmoidal & is.null(t)) stop('Provide non-null t for pseudotime model fitting')
#   model <- NULL
#   if(dist == 'nb') {
#     if(sigmoidal) {
#       model <- fit_alt_model(x, t, obj_func = nb_alt_obj_func, dist = dist)
#     } else {
#       model <- fit_null_model(x, obj_func = nb_null_obj_func, dist = dist)
#     }
#   } else if(dist == "normal") {
#     if(sigmoidal) {
#       model <- fit_alt_model(x, t, obj_func = norm_alt_obj_func, dist = dist)
#     } else {
#       model <- fit_null_model(x, obj_func = norm_null_obj_func, dist = dist)
#     }
#   }
#   model$type <- ifelse(sigmoidal, "sigmoid", "flat")
#   model$dist <- dist
#   return(model)
# }


lrtest <- function(x, t, pst_model, null_model) {
  ## first alternative model
  if(is.na(pst_model)) return(1)
  if(is.na(null_model)) return(0)
  pst_log_lik <- pst_model$log_lik
  null_log_lik <- null_model$log_lik
  D <- 2 * (pst_log_lik - null_log_lik)
  dof <- pst_model$df - null_model$df
  return( pchisq(D, dof, lower.tail = FALSE))
}

diffExprTest <- function(x, t, dist = c('nb', 'normal')) {
  dist <- match.arg(dist)
  pst_model <- fitModel(x, t, dist = dist, sigmoidal = TRUE)
  null_model <- fitModel(x, dist = dist, sigmoidal = FALSE)
  lrtest(x, t, pst_model, null_model)
}



plot_model <- function(model, x, t) {

  df <- data.frame(y=x, t=t)
  plt <- ggplot(df) + geom_point(aes(x=t, y=y), alpha=.8, size = 3) +
    xlab('Pseudotime') + ylab('Expression')

  if(model$type == "sigmoid") {
    mu_func <- function(t, params) {
      L <- params[1] ; k <- params[2] ; t_0 <- params[3] ; r <- params[4]
      L / (1 + exp(-k*(t - t_0)))
    }
    plt <- plt + stat_function(fun = mu_func, args = list(model$par), color='red', size = 1.5)
  } else if(model$type == "flat") {
    mu_func <- function(t, params) {
      params[1]
    }
    plt <- plt + stat_function(fun = mu_func, args = list(model$par), color='red', size = 1.5)
  }

  if ( library(cowplot, logical.return = TRUE) ) {
    plt <- plt + cowplot::theme_cowplot()
  } else {
    plt <- plt + theme_bw()
  }
  return( plt )

}

#' Differential expression test across pseudotime
diffExprSCE <- function(sce, use = "expr", dist = c('nn','normal')) {
  if(is.null(sce$pseudotime) && is.null(sce$Pseudotime)) {
    stop('Fit pseudotime before differential expression test')
  }
  X <- switch(use,
          expr = exprs(sce),
          tpm = tpm(sce),
          fpkm = fpkm(sce))

  pt <- pseudotime(sce)

  pvals <- apply(X, 1, diffExprTest, pt, dist)

  df_pval <- data.frame(pval = pvals)
  df_pval$gene <- row.names(df_pval)
  df_pval$qval <- p.adjust(df_pval$pval, method = 'BH')
  return(df_pval)
}

#' Fit pseudotime model hypers
#'
#' @importFrom MASS fitdistr
fitPstModelHypers <- function(model_list, trim = 0.05) {
  pst_params <- do.call("rbind", lapply(model_list, '[[', 'par'))
  pst_params_trimmed <- apply(pst_params, 2, function(p) {
    qnt <- quantile(p, probs = c(trim, 1 - trim))
    to_keep <- p > qnt[1] & p < qnt[2]
    p[to_keep]
  })
  ## going to fit gamma normal normal gamma
  L_fit <- fitdistr(pst_params_trimmed[,1], densfun = 'gamma')
  k_fit <- c(mu = mean(pst_params_trimmed[,2]), sd = sd(pst_params_trimmed[,2]))
  t0_fit <- c(mu = mean(pst_params_trimmed[,3]), sd = sd(pst_params_trimmed[,3]))
  r_fit <- fitdistr(pst_params_trimmed[,4], densfun = 'gamma')
  return(list(L = L_fit$estimate,
              k = k_fit,
              t0 = t0_fit,
              r = r_fit$estimate))
}

fitNullModelHypers <- function(model_list, trim = 0.05) {
  pst_params <- do.call("rbind", lapply(model_list, '[[', 'par'))
  pst_params_trimmed <- apply(pst_params, 2, function(p) {
    qnt <- quantile(p, probs = c(trim, 1 - trim))
    to_keep <- p > qnt[1] & p < qnt[2]
    p[to_keep]
  })
  ## going to fit gamma normal normal gamma
  mu_fit <- fitdistr(pst_params_trimmed[,1], densfun = 'gamma')
  r_fit <- fitdistr(pst_params_trimmed[,2], densfun = 'gamma')
  return(list(mu = mu_fit$estimate,
              r = r_fit$estimate))
}
