


#' @title Kolmogorov-Smirnov Tests for Distribution Estimates
#' 
#' @description 
#' Perform the Kolmogorov-Smirnov tests for various distribution estimates, via \code{\link[stats]{ks.test}}.
#' 
#' @param x an R object of distribution estimates
#' 
#' @param ... additional parameters of \code{\link[stats]{ks.test}}
#' 
#' @return 
#' 
#' \code{\link{ks_test}} returns the \code{'htest'} object (returned from \code{\link[stats]{ks.test}}),
#' in which the element \code{$statistic} is the Kolmogorov–Smirnov distance.
#' 
#' @seealso \code{\link[stats]{ks.test}}
#' 
#' @export
ks_test <- function(x, ...) UseMethod('ks_test')

#' @export
ks_test.fmx_QLMDe <- function(x, ...) {
  ks.test(x = x@data, y = pfmx, dist = x, ...)
}

#' @export
ks_test.fitdist <- function(x, ...) {
  if (!length(x$data)) stop('Re-run ?fitdistrplus::fitdist with option `keepdata = TRUE`')
  do.call(ks.test, args = c(list(x = quote(x$data), y = paste0('p', x$distname), ...), as.list.default(x$estimate)))
}


