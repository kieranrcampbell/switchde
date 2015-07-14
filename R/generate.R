## generate synthetic pseudotemporally regulated scRNA-seq data

generate_pst_cells <- function(ngenes, hypers, pst, dropout = TRUE, rate = NULL, base = 10) {
  L <- rgamma(ngenes, shape = hypers$L[1], rate = hypers$L[2])
  k <- rnorm(ngenes, mean = hypers$k[1], sd = hypers$k[2])
  t0 <- rnorm(ngenes, mean = hypers$t0[1], sd = hypers$t0[2])
  r <- rgamma(ngenes, shape = hypers$r[1], rate = hypers$r[2])
  param_matrix <- cbind(L, k, t0, r)
  mus <- apply(param_matrix, 1, calc_mu, pst)
  X <- sapply(1:ngenes, function(i) {
    y <- rnbinom(length(pst), mu = mus[,i], size = r[i])
    if(dropout) {
      x <- log(mus[,i] + 1, base = base)
      p_dropout <- exp(-rate * x * x)
      is_amplified <- sapply(p_dropout, function(p) sample(c(0,1), size = 1, prob = c(p, 1-p)))
      y <- y * is_amplified
    }
    y
  })
  return( X )
}

generate_null_cells <- function(ngenes, ncells, hypers, dropout = TRUE, rate = NULL, base = 10) {

}

#' Simulate single-cell RNA-seq data
#'
#' @export
simulateCells <- function(ngenes, hypers, pst = NULL,
                          ncells = NULL, dropout = TRUE, rate = NULL, base = 10) {
  if(dropout && is.null(rate)) stop('Specify rate for dropout = TRUE')
  if(is.null(pst) && is.null(ncells)) stop('If no pseudotime specify number of cells')
  X <- NULL
  if(!is.null(pst)) {
    X <- generate_pst_cells(ngenes, hypers, pst, dropout, rate)
  } else {
    X <-  generate_pst_cells(ngenes, ncells, hypers, dropout, rate)
  }
  if(is.null(ncells)) ncells <- length(pst)
  colnames(X) <- paste0('gene', 1:ngenes)
  rownames(X) <- paste0('cell', 1:ncells)
  return(X)
}
