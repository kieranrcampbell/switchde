

#' Switch-like model fitting and differential expression test
#' 
#' Fit sigmoidal differential expression models to gene expression across pseudotime.
#' Parameter estimates are returned along with a p-value for switch-like differential
#' expression over a null model (constant expression).
#'  
#' @param object Gene expression data that is either
#' \itemize{
#'  \item A vector of length number of cells for a single gene
#'  \item A matrix of dimension number of genes x number of cells
#'  \item An object of class \code{SingleCellExperiment} from package SingleCellExperiment
#'  }
#' @param pseudotime A pseudotime vector with a pseudotime corresponding to 
#' every cell. Can be \code{NULL} if object is of class \code{SCESet} and 
#' \code{pData(sce)$pseudotime} is defined.
#' @param zero_inflated Logical. Should zero inflation be implemented? Default  \code{FALSE}
#' @param lower_threshold The minimum threshold below which to set expression to zero to avoid
#' numerical issues. Default is 0.01
#' @param maxiter Maximum number of iterations for EM algorithm if zero inflation enabled. Default 100
#' @param log_lik_tol If the change in the log-likelihood falls below this for zero inflated EM
#' the algorithm is assumed to have converged
#' @param verbose Print convergence update for EM algorithm
#' @param sce_assay The assay from the \code{SingleCellExperiment} to be used
#' as expression, defaulting to "exprs"
#'  
#' @export
#' 
#' @import dplyr
#' 
#' @return A matrix where each column corresponds to a gene, the first row is
#' the p-value for that gene and subsequent rows are model parameters.
#' 
#' @importFrom stats p.adjust
#' 
#' @examples
#' data(synth_gex)
#' data(ex_pseudotime)
#' sde <- switchde(synth_gex, ex_pseudotime)
switchde <- function(object, pseudotime = NULL, zero_inflated = FALSE,
                     lower_threshold = 0.01, maxiter = 1000, 
                     log_lik_tol = 1e-2, verbose = FALSE,
                     sce_assay = "exprs") {
  res <- NULL
  inputs <- sanitise_inputs(object, pseudotime, lower_threshold, zero_inflated, sce_assay)
  X <- inputs$X
  pst <- inputs$pst
  
  feature_names <- NULL
  if(!is.null(rownames(X))) {
    feature_names <- rownames(X)
  } else {
    feature_names <- paste0("gene", 1:nrow(X))
  }
  
  if(any(is.na(inputs$pst))) {
    warning("NAs present in pseudotime vector. Returning NA for all parameters.")
    na_vec <- rep(NA, nrow(X))
    res <- dplyr::data_frame(gene = feature_names, 
                             pval = na_vec, qval = na_vec, 
                             mu0 = na_vec, k = na_vec, t0 = na_vec)
    if(zero_inflated) {
      res <- mutate(res, lambda = na_vec, EM_converged = na_vec)
    }
    
    return(res)
  }
  
  
  ## differential gene test time
  if(zero_inflated) {
    res <- apply(X, 1, fit_zi_model, pst, maxiter, log_lik_tol, verbose)
  } else {
    res <- apply(X, 1, fit_nzi_model, pst)
  }
  
  res <- as_data_frame(t(res))
  res <- mutate(res, gene = feature_names)
  
  
  ## This just appeases R CMD CHECK
  pval <- gene <- qval <- mu0 <- k <- t0 <- lambda <- EM_converged <- NULL
  
  res <- mutate(res, qval = p.adjust(pval, method = "BH"))
  
  if(zero_inflated) {
    res <- select(res, gene, pval, qval, mu0, k, t0, lambda, EM_converged)
    res$EM_converged <- as.logical(res$EM_converged)
    if(any(!res$EM_converged)) {
      warning("Some EM fits did not converge. Please inspect the `EM_converged` column")
    }
  } else {
    res <- select(res, gene, pval, qval, mu0, k, t0)
  }
  
  return( res )
}

#' Extract parameters from fitted model
#' 
#' Extract maximum likelihood parameter estimates from a call to \code{switchde}.
#' 
#' @param sde The \code{data.frame} returned by \code{switchde}
#' @param gene The gene for which to extract parameters
#' @return A vector of length 3 corresonding to the parameters \eqn{\mu_0}, \eqn{k} and \eqn{t_0}
#' 
#' 
#' @export
#' @examples
#' data(synth_gex)
#' data(ex_pseudotime)
#' sde <- switchde(synth_gex, ex_pseudotime)
#' pars <- extract_pars(sde, "Gene1")
extract_pars <- function(sde, gene) {
  stopifnot(gene %in% sde$gene)
  g <- gene
  sde_gene <- filter(sde, gene == g) 
  return( unlist(sde_gene[,4:6]) )
}

#' Plot gene behaviour
#' 
#' Plot gene behaviour and MLE sigmoid as a function of pseudotime.
#' 
#' @param x Gene expression vector
#' @param pseudotime Pseudotime vector (of same length as x)
#' @param pars Fitted model parameters
#' 
#' @details This plots expression of a single gene. Fitted model parameters can
#' either be specified manually or can be extracted from the \code{data.frame} returned
#' by \code{switchde} using the function \code{extract_pars}.
#' 
#' @import ggplot2
#' @export
#' 
#' @return A \code{ggplot2} plot of gene expression and MLE sigmoid
#' 
#' @examples
#' data(synth_gex)
#' data(ex_pseudotime)
#' sde <- switchde(synth_gex, ex_pseudotime)
#' switchplot(synth_gex[1, ], ex_pseudotime, extract_pars(sde, "Gene1"))
switchplot <- function(x, pseudotime, pars) {
  ggplot(data_frame(Expression = x, Pseudotime = pseudotime), 
         aes_string(x = "Pseudotime", y = "Expression")) +
    geom_point(alpha = 0.5, fill = "grey", colour = "black", shape = 21) + theme_bw() +
    stat_function(fun = sigmoid, args = list(params = pars), color = 'red')
}



#' Sanitise inputs 
#' @param object The object passed at the entry point (either a SCESet or gene
#' expression matrix)
#' @param pseudotime A pseudotime vector
#' @param zero_inflated Logical. Should zero inflation be implemented? Default  \code{FALSE}
#' @param lower_threshold The minimum threshold below which to set expression to zero to avoid
#' numerical issues. Default is 0.01
#' @param sce_assay The assay from the \code{SingleCellExperiment} to be used
#' as expression, defaulting to "exprs"
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @keywords internal
#' 
#' @return A list with two entries: a gene expression matrix \code{X}
#' and a pseudotime vector \code{pst}.
sanitise_inputs <- function(object, pseudotime, lower_threshold, zero_inflated,
                            sce_assay) {
  X <- pst <- NULL
  
  if(is.vector(object) && is.numeric(object)) { # single gene expression vector
    message(paste("Assuming single gene measured in", length(object)), " cells")
    X <- matrix(object, nrow = 1)
  } else if(is.matrix(object)) { # multiple gene expression matrix
    message(paste("Input gene-by-cell matrix:", nrow(object), "genes and ", ncol(object), "cells"))
    X <- object
  } else if(is(object, "SingleCellExperiment")) {
    if(is.null(pseudotime)) {
      pst <- pData(object)$pseudotime
    }
    if(is.null(pseudotime)) stop("Pseudotime must either be specified or as a column named pseudotime in the phenoData of the SCESet")
    X <- assay(object, sce_assay)
  } else {
    stop("object must be vector, matrix or SingleCellExperiment")
  }
  
  if(is.null(pst) && !is.null(pseudotime)) pst <- pseudotime
  if(is.null(X)) stop("Object must either be numeric vector, matrix or SingleCellExperiment")
  if(is.null(pst)) stop("Pseudotime must either be specified or in pData(object)")
  if(length(pst) != ncol(X)) stop("Must have pseudotime for each cell")
  
  ## lower threshold - some calculations in EM suffer numerical issues for
  ## very small y which are meaningless anyway
  X[X < lower_threshold] <- 0
  
  gene_vars <- apply(X, 1, var)
  num_genes_zero_var <- sum(gene_vars == 0)
  if(num_genes_zero_var > 0) {
    msg <- paste("After setting measurements below lower_threshold =", lower_threshold)
    msg <- paste(msg, "to zero,", num_genes_zero_var, "features had zero variance.")
    msg <- paste(msg, "Please consider stricter filtering, such as mean expression above")
    msg <- paste(msg, "lower_threshold and nonzero expression in 10% of cells")
    stop(msg)
  }
  
  
  
  return( list(X = X, pst = pst) )
}


#' Example sigmoid plot
#' 
#' Plot an example sigmoid function. For demonstration and documentation.
#' 
#' 
#' @import ggplot2
#' @export
#' 
#' @return An object of class \code{ggplot}
#' 
#' @examples
#' example_sigmoid()
example_sigmoid <- function() {
  
  theme_set(theme_bw())
  
  f <- function(x, k, mu0, t0) {
    2 * mu0 / (1 + exp(-k * (x - t0)))
  }
  
  g <- function(x, m, C) {
    m * x + C
  }
  
  k <- 20
  mu0 <- 1
  t0 <- 0.5
  
  cols = c("#E41A1C", "#377EB8", "#4DAF4A")
  
  plt <- ggplot(data.frame(x = c(0,1)), aes_string(x = "x")) +
    theme_bw() + 
    xlab("Pseudotime")
  
  plt <- plt + geom_hline(yintercept = mu0, colour = cols[2], linetype = 1, alpha = 0.8, size = 1.5) +
    geom_vline(xintercept = 0.5, colour = cols[3], linetype = 1, size = 1.5, alpha = 0.8) +
    stat_function(fun = g, args = list(C  = -3.5, m = 9), 
                  linetype = 1, alpha = 0.8, colour = cols[1], size = 1.5)
  
  plt <- plt + geom_text(label = "mu[0]", x = 0.25, y = 1.3, parse = TRUE,
                         size = 8, colour = cols[2]) +
    geom_text(label = "t[0]", x = 0.55, y = 0, parse = TRUE,
              size = 8, colour = cols[3]) +
    geom_text(label = "k", x = 0.55, y = 2, size = 8, colour = cols[1])
  
  plt <- plt + stat_function(fun = f, args = list(k = k, mu0 = mu0, t0 = t0), size = 1.5, alpha = 0.7) +
    ylim(-0.5,2.5) + ylab("Gene expression")
  
  return(plt)
}

#' Fit a zero-inflated model for a single gene
#' 
#' Fits a zero-inflated sigmoidal model for a single gene vector, returning
#' MLE model parameters and p-value.
#' 
#' @param y Vector of gene expression values
#' @param pst Pseudotime vector, of same length as y
#' @param maxiter Maximum number of iterations for EM algorithm if zero inflation enabled. Default 100
#' @param log_lik_tol If the change in the log-likelihood falls below this for zero inflated EM
#' the algorithm is assumed to have converged
#' @param verbose Print convergence update for EM algorithm
#' 
#' @export
#' @return A vector with 6 entries: maximum likelihood estimates for \eqn{\mu_0}, \eqn{k}
#' \eqn{t0}, \eqn{\lambda}, \eqn{\sigma^2} and a p-value
#' 
#' @importFrom stats pchisq
#' 
#' @examples
#' data(synth_gex)
#' data(ex_pseudotime)
#' y <- synth_gex[1, ]
#' fit <- fit_zi_model(y, ex_pseudotime)
fit_zi_model <- function(y, pst, maxiter = 10000, log_lik_tol = 1e-3, verbose = FALSE) {
  
  stopifnot(length(y) == length(pst))
  stopifnot(all(y >= 0))
  
  if(var(y) == 0) stop("Variance of input expression is zero - cannot fit temporal model.")
  
  if(!any(y == 0)) { # Return standard model
    nzi_results <- fit_nzi_model(y, pst)
    r <- rep(NA, 7)
    names(r) <- c("mu0", "k", "t0", "sigma2", "lambda", "pval", "EM_converged")
    r[c(1:4,6)] <- nzi_results
    r[7] <- TRUE
    return(r)
  }
  
  sigmoidal_model <- EM_sigmoid(y, pst, iter = maxiter, 
                                log_lik_tol = log_lik_tol, verbose = verbose)
  
  
  r <- NULL
  
  constant_model <- EM_constant(y, maxiter, log_lik_tol, verbose)
  
  D <- -2 * (constant_model$log_lik - sigmoidal_model$log_lik)
  dof <- 2 # 4 - 2
  pval <- pchisq(D, dof, lower.tail = FALSE)
  
  r <- c(sigmoidal_model$params, pval, 
         (1 * sigmoidal_model$converged) * constant_model$converged)
  
  
  names(r) <- c("mu0", "k", "t0", "sigma2", "lambda", "pval", "EM_converged")
  return(r)
}

#' Fit a (non-zero-inflated) model for a single gene
#' 
#' Fits a sigmoidal expression model for a single gene vector, returning
#' MLE model parameters and p-value.
#' 
#' @param y Vector of gene expression values
#' @param pst Pseudotime vector, of same length as y
#' 
#' @export
#' @return A vector with 5 entries: maximum likelihood estimates for \eqn{\mu_0}, \eqn{k}
#' \eqn{t0}, \eqn{\sigma^2} and a p-value
#' 
#' @examples
#' data(synth_gex)
#' data(ex_pseudotime)
#' y <- synth_gex[1, ]
#' fit <- fit_nzi_model(y, ex_pseudotime)
fit_nzi_model <- function(y, pst) {
  sigmoidal_model <- fit_sigmoidal_model(y, pst)
  constant_model <- fit_constant_model(y)
  
  D <- -2 * (constant_model$log_lik - sigmoidal_model$log_lik)
  dof <- 2 # 4 - 2
  pval <- pchisq(D, dof, lower.tail = FALSE)
  
  r <- c(sigmoidal_model$params, pval)
  names(r) <- c("mu0", "k", "t0", "sigma2", "pval")
  return(r)
}

#' Calculate the mean vector given parameters and pseudotimes (mu0 formulation)
#' 
#' This function (common to all models) calculates the sigmoidal mean vector
#' given the parameters and factor of pseudotimes
#' 
#' @param params Vector of length 3 with entries mu_0, k, t0
#' @param pst Vector of pseudotimes
#' 
#' @keywords internal
#' 
#' @return Mean sigmoidal vector
sigmoid <- function(pst, params) {
  mu0 <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- 2 * mu0 / (1 + exp(-k*(pst - t_0)))
  return(mu)
}


#' Synthetic gene expression matrix
#' 
#' A matrix containing some synthetic gene expression data for 
#' 100 cells and 12 genes
#' 
#' @return A 12 by 100 matrix
"synth_gex"

#' Synthetic gene pseudotimes
#' 
#' A vector with example pseudotimes for the synthetic 
#' gene expression data in \code{example_gex}
#' @return A vector of length 100
"ex_pseudotime"