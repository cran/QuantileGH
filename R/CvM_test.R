


#' @title Cramer-Von Mises Test of Goodness-of-Fit for Distribution Estimates
#' 
#' @description
#' 
#' Perform the Cramer-Von Mises test of goodness-of-fit for distribution estimates, via \code{\link[goftest]{cvm.test}}.
#' 
#' @param x an R object of distribution estimates
#' 
#' @param nullname,... additional parameters of \code{\link[goftest]{cvm.test}}
#' 
#' @details 
#' 
#' Note that we are currently not using the discrete version \code{\link[dgof]{cvm.test}}.
#' 
#' @return \code{\link[goftest]{cvm.test}}
#' 
#' \code{\link{CvM_test}} returns the \code{'htest'} object (returned from \code{\link[goftest]{cvm.test}}),
#' in which the element \code{$statistic} is the Cramer-Von Mises quadratic distance.
#' 
#' @seealso \code{\link[goftest]{cvm.test}}
#' 
#' @export
CvM_test <- function(x, nullname, ...) UseMethod('CvM_test')

#' @export
CvM_test.fmx_QLMDe <- function(x, nullname = deparse1(substitute(x)), ...) {
  # cvm.test(x = x@data, null = pfmx, dist = x, ...) # causes error in \strong{printing} !!
  cvm.test(x = x@data, null = function(q) pfmx(q, dist = x), nullname = nullname, ...)
}

#' @export
CvM_test.fitdist <- function(x, nullname = deparse1(substitute(x)), ...) {
  if (!length(x$data)) stop('Re-run ?fitdistrplus::fitdist with option `keepdata = TRUE`')
  cvm.test(x = x$data, null = function(q) do.call(paste0('p', x$distname), args = c(list(q = q), as.list.default(x$estimate))), nullname = nullname, ...)
}





