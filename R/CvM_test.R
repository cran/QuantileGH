


#' @title Cramer-Von Mises Test of Goodness-of-Fit for Distribution Estimates
#' 
#' @description
#' 
#' Perform the Cramer-Von Mises test of goodness-of-fit for distribution estimates.
#' 
#' @param x an R object of distribution estimates
#' 
#' @param data \link[base]{double} vector, the actual observations
#' 
#' @param nullname,... additional parameters of \code{\link[goftest]{cvm.test}}
#' 
#' @details 
#' 
#' Note that we are currently not using the discrete version \code{\link[dgof]{cvm.test}}.
#' 
#' @return 
#' 
#' \link{CvM_test} returns an \link[goftest:cvm.test]{htest} object,
#' in which the element \code{$statistic} is the Cramer-Von Mises quadratic distance.
#' 
#' @seealso \link[goftest]{cvm.test}
#' 
#' @examples 
#' (d1 <- fmx('norm', mean = c(0, 1.5), sd = .5, w = c(.4, .6)))
#' x = rfmx(1e2L, dist = d1)
#' CvM_test(d1, data = x)
#' 
#' @export
CvM_test <- function(x, data, nullname, ...) UseMethod('CvM_test')

#' @export
CvM_test.fmx <- function(x, data = x@data, nullname = deparse1(substitute(x)), ...) {
  if (!length(data)) stop('must provide actual observations')
  # cvm.test(x = data, null = pfmx, dist = x, nullname = nullname, ...) # cannot deal complex parameters such as my 'fmx'
  cvm.test(x = data, null = function(q) {
    pmin.int(pmax.int(pfmx(q, dist = x), .Machine$double.eps), 1 - .Machine$double.eps)
  }, nullname = nullname, ...)
}

#' @export
CvM_test.fitdist <- function(x, data = x$data, nullname = deparse1(substitute(x)), ...) {
  if (!length(data)) stop('Re-run ?fitdistrplus::fitdist with option `keepdata = TRUE`')
  cvm.test(x = data, null = function(q) {
    tmp <- do.call(paste0('p', x$distname), args = c(list(q = q), as.list.default(x$estimate)))
    pmin.int(pmax.int(tmp, .Machine$double.eps), 1 - .Machine$double.eps)
  }, nullname = nullname, ...)
}





