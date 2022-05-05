


#' @title Kolmogorov-Smirnov Tests for Distribution Estimates
#' 
#' @description 
#' Perform the Kolmogorov-Smirnov tests for various distribution estimates.
#' 
#' @param x an R object of distribution estimates
#' 
#' @param data \link[base]{double} vector, the actual observations
#' 
#' @param ... additional parameters of \link[stats]{ks.test}
#' 
#' @return 
#' 
#' \link{ks_test} returns an \link[stats:ks.test]{htest} object,
#' in which the element \code{$statistic} is the Kolmogorov–Smirnov distance.
#' 
#' @seealso \link[stats]{ks.test}
#' 
#' @export
ks_test <- function(x, data, ...) UseMethod('ks_test')

#' @export
ks_test.fmx <- function(x, data = x@data, ...) {
  if (!length(data)) stop('must provide actual observations')
  ks.test(x = data, y = pfmx, dist = x, ...)
}

#' @export
ks_test.fitdist <- function(x, data = x$data, ...) {
  if (!length(data)) stop('Re-run ?fitdistrplus::fitdist with option `keepdata = TRUE`')
  do.call(ks.test, args = c(list(x = quote(data), y = paste0('p', x$distname), ...), as.list.default(x$estimate)))
}


