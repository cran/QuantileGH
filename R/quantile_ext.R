


#' @title Variance-Covariance of Quantiles
#' 
#' @description
#' 
#' Computes the variance-covariance matrix of quantiles based on Theorem 1 and 2 of Mosteller (1946).
#' 
#' @param probs \link[base]{numeric} \link[base]{vector}, cumulative probabilities at the given quantiles
#' 
#' @param d \link[base]{numeric} \link[base]{vector}, probability densities at the given quantiles
#' 
#' @details 
#' 
#' The end user should make sure no density too close to 0 is included in argument `d`.
#' 
#' Function [quantile_vcov()] must not be used in a compute-intensive way.
#' 
#' @returns 
#' 
#' Function [quantile_vcov()] returns the variance-covariance \link[base]{matrix} of quantiles based on Mosteller (1946).
#' 
#' @references 
#' Frederick Mosteller. On Some Useful "Inefficient" Statistics (1946).
#' \doi{10.1214/aoms/1177730881}
#' 
#' @export
quantile_vcov <- function(probs, d) {
  # do the check on d=0 in \link{QLMDe}
  if (anyNA(probs) || anyNA(d)) stop('no NA allowed in probability nor density')
  if ((n <- length(probs)) != length(d)) stop('probs and d must match in length')
  
  fs <- tcrossprod(d, d) # 'matrix'
  p_c <- array(probs, dim = c(n,n)) # 'probs on cols'
  p_r <- t.default(p_c) # 'probs on rows'
  p_min <- pmin.int(p_r, p_c) # vector!
  p_max <- pmax.int(p_r, p_c)
  vv <- p_min * (1 - p_max) / fs # back to 'matrix'
  return(vv)
}