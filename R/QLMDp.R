
#' @title Percentages for Quantile Least Mahalanobis Distance estimation
#' 
#' @description 
#' 
#' A vector of probabilities to be used in Quantile Least Mahalanobis Distance estimation (\code{\link{QLMDe}}). 
#' 
#' @param from \link[base]{numeric} scalar, minimum of the equidistant (in probability or quantile) probabilities.  Default \code{.05}.
#' 
#' @param to \link[base]{numeric} scalar, maximum of the equidistant (in probability or quantile) probabilities.  Default \code{.95}.
#' 
#' @param length.out non-negative \link[base]{integer} scalar, the number of the equidistant (in probability or quantile) probabilities. 
#' 
#' @param equidistant \link[base]{character} scalar.
#' If \code{'prob'} (default), then the probabilities are equidistant.  
#' If \code{'quantile'}, then the quantiles (of the observations \code{x}) corresponding to the probabilities are equidistant.
#' 
#' @param extra \link[base]{numeric} vector of \emph{additional} probabilities, default \code{c(.005, .01, .02, .03, .97, .98, .99, .995)}.
#' 
#' @param x \link[base]{numeric} vector of observations, only used when \code{equidistant = 'quantile'}.
#' 
#' @details
#' 
#' The default arguments of \link{QLMDp} returns the probabilities of 
#' \code{c(.005, .01, .02, .03, seq.int(.05, .95, length.out = 15L), .97, .98, .99, .995)}.
#' 
#' @return 
#' 
#' A \link[base]{numeric} vector of probabilities to be supplied to parameter \code{p} of 
#' Quantile Least Mahalanobis Distance \link{QLMDe} estimation).
#' In practice, the length of this probability vector \code{p} 
#' must be equal or larger than the number of parameters in the distribution model to be estimated.
#' 
#' @examples 
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(100); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' p_hist = geom_histogram(
#'   mapping = aes(x = x2, y = ..density..), bins = 30L, colour = 'white', alpha = .1
#' )
#' 
#' (p1 = QLMDp()) # equidistant in probabilities
#' autoplot(d2, v = setNames(qfmx(p1, dist = d2), nm = sprintf('%.1f%%', 1e2*p1)))
#' autoplot(d2, v = quantile(x2, probs = p1, digits = 3L)) + p_hist
#' 
#' (p2 = QLMDp(equidistant = 'quantile', x = x2)) # equidistnat in quantiles
#' autoplot(d2, v = quantile(x2, probs = p2, digits = 3L)) + p_hist
#' 
#' 
#' @export
QLMDp <- function(
  from = .05, to = .95, length.out = 15L, equidistant = c('prob', 'quantile'),
  extra = c(.005, .01, .02, .03, .97, .98, .99, .995),
  x
) {
  
  if (!is.double(from) || length(from) != 1L || is.na(from) || from <= 0) stop('`from` must be numeric >0')
  if (!is.double(to) || length(to) != 1L || is.na(to) || to >= 1) stop('`to` must be numeric <1')
  if (!is.integer(length.out) || length(length.out) != 1L || anyNA(length.out) || length.out < 0L) stop('N must be len-1 non-negative integer')
  if (!is.double(extra) || anyNA(extra) || any(extra <= 0, extra >= 1)) stop('`extra` must be numeric between (0, 1)') # complatible with len-0

  p <- if (!length.out) extra else {
    c(extra, switch(match.arg(equidistant), prob = {
      #cat('Equidistant in probability; ')
      seq.int(from = from, to = to, length.out = length.out)
    }, quantile = {
      #cat('Equidistant in quantile; ')
      if (missing(x)) stop('Must have `x`')
      qlim <- quantile(x, probs = c(from, to)) # not compute intensive
      ecdf(x)(seq.int(from = qlim[1L], to = qlim[2L], length.out = length.out))
    }))
  }
  if (!length(p)) stop('len-0 p, do not allow')
  sort.int(unique_allequal(p))
}




