


#' @title Kullback-Leibler Divergence for Distribution Estimates
#' 
#' @description 
#' 
#' Calculate the Kullback-Leibler divergence for distribution estimates, via \code{\link[LaplacesDemon]{KLD}}.
#' 
#' @param x an R object of distribution estimates
#' 
#' @param base see \code{\link[LaplacesDemon]{KLD}}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @return 
#' 
#' \code{\link{KL_dist}} returns a \code{\link[base]{list}}, which is returned from \code{\link[LaplacesDemon]{KLD}}.
#' 
#' @seealso \code{\link[LaplacesDemon]{KLD}}.
#' 
#' @export
KL_dist <- function(x, base, ...) UseMethod('KL_dist')

#' @export
KL_dist.fmx_QLMDe <- function(x, base = exp(1), ...) {
  KLD(px = x@epdf(x@data), py = dfmx(x@data, dist = x, log = FALSE), base = base)
}

#' @export
KL_dist.fitdist <- function(x, base = exp(1), ...) {
  if (!length(x$data)) stop('Re-run ?fitdistrplus::fitdist with option `keepdata = TRUE`')
  x_kern <- density.default(x$data) # ?stats::approx inside ?stats::density.default
  x_epdf <- approxfun(x = x_kern$x, y = x_kern$y) # another 'layer' of ?stats::approxfun
  KLD(px = x_epdf(x$data), py = do.call(paste0('d', x$distname), args = c(list(x = x$data, log = FALSE), as.list.default(x$estimate))), base = base)
}