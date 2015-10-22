# Sigmoidal differential expression for single-cell
# RNA-seq data across pseudotime.
# 
# kieran.campbell@sjc.ox.ac.uk

#' Returns a set of fitted sigmoidal pseudotime models
#' 
#' @export
fitModel <- function(object, ...) {
  
}

#' Tests differential expression for one or more genes
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
#' @return A matrix where each column corresponds to a gene, the first row is
#' the p-value for that gene and subsequent rows are model parameters.
testDE <- function(object, pseudotime = NULL, zero_inflated = FALSE, ...) {
  X <- pst <- res <- NULL
  
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
  
  ## differential gene test time
  if(zero_inflated) {
    res <- apply(X, 1, sctools::norm_zi_diff_expr_test, pst, ...)
  } else {
    res <- apply(X, 1, sctools::norm_diff_expr_test, pst)
  }
  
  return( res )
}




