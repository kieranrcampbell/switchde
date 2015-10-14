
# Common functions to all models fitted

#' Calculate the mean vector given parameters and pseudotimes
#' 
#' This function (common to all models) calculates the sigmoidal mean vector
#' given the parameters and factor of pseudotimes
#' 
#' @param params Vector of length 4 with entries L = 2\mu_0, k, t0, sig_sq
#' @param t Vector of pseudotimes
#' 
#' @return Mean sigmoidal vector
#' 
#' @export
calc_mu <- function(params, t) {
  L <- params[1] ; k <- params[2] ; t_0 <- params[3]
  mu <- L / (1 + exp(-k*(t - t_0)))
  return(mu)
}
