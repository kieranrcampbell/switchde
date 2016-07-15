# Sigmoidal differential expression for single-cell
# RNA-seq data across pseudotime.
# 
# kieran.campbell@sjc.ox.ac.uk

#' Fit sigmoidal pseudotime model for set of genes
#' 
#'  @param object Gene expression data that is either
#'  \itemize{
#'  \item A vector of length number of cells for a single gene
#'  \item A matrix of dimension number of genes x number of cells
#'  \item An object of class \code{SCESet} from package scater
#'  }
#'  @param pseudotime A pseudotime vector with a pseudotime corresponding to 
#'  every cell. Can be \code{NULL} if object is of class \code{SCESet} and 
#'  \code{pData(sce)$pseudotime} is defined.
#'  @param zero_inflated Logical. Should zero inflation be implemented? Default  \code{FALSE}
#'  @param ... Additional arguments to be passed to expectation maximisation algorithm
#'  if zero-inflation is enabled:
#'  \itemize{
#'  \item maxit Maximum number of iterations for EM. Default 500
#'  \item loglik_tol Change in log-likelihood tolerance. Default 1e-7 
#'  \item verbose Should progress of EM optimisation be printed? Default \code{FALSE}
#'  }
#'  
#' @export
#' 
#' @import dplyr
#' 
#' @return A matrix where columns are genes and rows are model parameters (zero inflation will
#' have an extra lambda parameter compared to no zero-inflation)
fitModel <- function(object, pseudotime = NULL, zero_inflated = FALSE, ...) {
  res <- NULL
  inputs <- sanitise_inputs(object, pseudotime)
  X <- inputs$X
  pst <- inputs$pst
  
  if(zero_inflated) {
    res <- apply(X, 1, function(x, pst, ...) {
      alt_EM_ctrl(x, pst, ...)$par
      }, pst, ...)
  } else {
    res <- apply(X, 1, function(x, pst) {
      norm_fit_alt_model(x, pst)$par
      }, pst)
  }
  
  return( res )
}

#' Switch-like differential expression test
#'  
#' @param object Gene expression data that is either
#' \itemize{
#'  \item A vector of length number of cells for a single gene
#'  \item A matrix of dimension number of genes x number of cells
#'  \item An object of class \code{SCESet} from package scater
#'  }
#' @param pseudotime A pseudotime vector with a pseudotime corresponding to 
#' every cell. Can be \code{NULL} if object is of class \code{SCESet} and 
#' \code{pData(sce)$pseudotime} is defined.
#' @param zero_inflated Logical. Should zero inflation be implemented? Default  \code{FALSE}
#' @param ... Additional arguments to be passed to expectation maximisation algorithm
#' if zero-inflation is enabled:
#' \itemize{
#'  \item maxit Maximum number of iterations for EM. Default 500
#'  \item loglik_tol Change in log-likelihood tolerance. Default 1e-7 
#'  \item verbose Should progress of EM optimisation be printed? Default \code{FALSE}
#'  }
#'  
#' @export
#' 
#' @return A matrix where each column corresponds to a gene, the first row is
#' the p-value for that gene and subsequent rows are model parameters.
switchde <- function(object, pseudotime = NULL, zero_inflated = FALSE, ...) {
  res <- NULL
  inputs <- sanitise_inputs(object, pseudotime)
  X <- inputs$X
  pst <- inputs$pst
  
  ## differential gene test time
  if(zero_inflated) {
    res <- apply(X, 1, norm_zi_diff_expr_test, pst, ...)
  } else {
    res <- apply(X, 1, norm_diff_expr_test, pst)
  }
  
  res <- tbl_df(t(res))
  
  if(!is.null(rownames(X))) res <- mutate(res, gene = rownames(X))
  
  res <- mutate(res, qval = p.adjust(pval, method = "BH"))
  res <- select(res, gene, pval, qval, mu0, k, t0)
  
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
#' @importFrom dplyr filter
#' 
#' @export
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
#' @param pst Pseudotime vector (of same length as x)
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
switchplot <- function(x, pseudotime, pars) {
  ggplot(data_frame(Expression = x, Pseudotime = pseudotime), aes(x = Pseudotime, y = Expression)) +
    geom_point(alpha = 0.5, fill = "grey", colour = "black", shape = 21) + theme_bw() +
    stat_function(fun = sigmoid, args = list(params = pars), color = 'red')
}

#' Backwards compatibility.
#' 
#' This function as been replaced by \code{switchde}.
#' 
#' @param ... Arguments passed to \code{switchde}.
#' 
#' @export
testDE <- function(...) {
  switchde(...)
}

#' Sanitise inputs for testDE and fitModel
sanitise_inputs <- function(object, pseudotime) {
  X <- pst <- NULL
  
  if(is.vector(object) && is.numeric(object)) { # single gene expression vector
    message(paste("Assuming single gene measured in", length(object)), "cells")
    X <- matrix(object, nrow = 1)
  } else if(is.matrix(object)) { # multiple gene expression matrix
    message(paste("Assuming gene-by-cell matrix:", nrow(object), "genes and ", ncol(object), "cells"))
    X <- object
  } else {
    if(require(scater)) {
      if(is(object, "SCESet")) {
        if(is.null(pseudotime)) pst <- pData(object)$pseudotime
        X <- exprs(object)
      } else {
        stop("object must be vector, matrix or SCESet")  
      }
    } else {
      stop("object must be vector, matrix or SCESet")
    }   
  }
  
  if(is.null(pst) && !is.null(pseudotime)) pst <- pseudotime
  if(is.null(X)) stop("Object must either be numeric vector, matrix or SCESet")
  if(is.null(pst)) stop("Pseudotime must either be specified or in pData(object)")
  if(length(pst) != ncol(X)) stop("Must have pseudotime for each cell")
  
  return( list(X = X, pst = pst) )
}


#' Example sigmoid plot
#' 
#' Mainly for documentation
#' 
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
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
  
  cols = brewer.pal(3, "Set1")
  
  plt <- ggplot(data.frame(x = c(0,1)), aes(x = x)) +
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