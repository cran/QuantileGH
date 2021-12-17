
#' @title Percentages used in Quantile Least Mahalanobis Distance estimation
#' 
#' @description 
#' 
#' The vector of probabilities corresponding to sample and theoretical quantiles selected for 
#' minimizing the Mahalanobis distance (\code{\link{QLMDe}}) is determined by function \code{\link{QLMDp}}.
#' 
#' @param p.lim range of the equidistant (in probability or quantile) probabilities.  Default \code{c(.05, .95)}.
#' 
#' @param p.N non-negative 'integer' value, the number of the equidistant (in probability or quantile) probabilities. 
#' \code{p.N = 0L} indicates no equidistant probabilities are to be used.
#' 
#' @param equidistant 'character' value, the type of equidistant probabilities to be used.  
#' If \code{'prob'} (default), then the probabilities are equidistant.  
#' If \code{'quantile'}, then the quantiles (of the observations \code{obs}) corresponding to the probabilities are equidistant.
#' 
#' @param p.extra \strong{additional} probabilities to be used, default \code{c(.005, .01, .02, .03, .05, .95, .97, .98, .99, .995)}.
#' 
#' @param obs 'numeric' vector of observations, only used when \code{equidistant = 'quantile'}
#' 
#' @param ... additional parameters
#' 
#' @details ..
#' 
#' @return A vector of probabilities, corresponding to which the sample and theoretical quantiles are selected for 
#' minimizing the Mahalanobis distance (\code{\link{QLMDe}}). 
#' In practice, the length of this probability vector must be equal or larger than the number of parameters in the distribution model
#' to be estimated by \code{\link{QLMDe}}.
#' 
#' @examples 
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(100); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' p_hist = geom_histogram(
#'   mapping = aes(x = x2, y = ..density..), 
#'   bins = 30L, colour = 'white', alpha = .1
#' )
#'   
#' (p1 = QLMDp()) # equidistant in probabilities
#' autoplot(d2, v = setNames(quantile(x2, probs = p1), nm = sprintf('%.1f%%', 1e2*p1))) + p_hist
#' 
#' (p2 = QLMDp(equidistant = 'quantile', obs = x2)) # equidistnat in quantiles
#' autoplot(d2, v = setNames(quantile(x2, probs = p2), nm = sprintf('%.1f%%', 1e2*p2))) + p_hist
#' 
#' 
#' @export
QLMDp <- function(
  p.lim = c(.05, .95),
  p.N = 15L,
  equidistant = c('prob', 'quantile'),
  p.extra = c(.005, .01, .02, .03, .05, .95, .97, .98, .99, .995),
  obs, # 'numeric' vector; observations 
  ...
) {
  
  if (!is.double(p.extra) || anyNA(p.extra) || any(p.extra <= 0, p.extra >= 1)) stop('`p.extra` must be numeric between (0, 1)') # complatible with len-0

  if (!is.double(p.lim) || anyNA(p.lim) || any(p.lim <= 0, p.lim >= 1)) stop('`p.lim` must be numeric between (0, 1)')
  if (length(p.lim) < 2L) stop('`p.lim` must be len >2L')
  p.lim <- c(min(p.lim), max(p.lim))
  
  if (!is.integer(p.N) || length(p.N) != 1L || anyNA(p.N) || p.N < 0L) stop('N must be len-1 non-negative integer')

  p <- if (!p.N) p.extra else {
    c(p.extra, switch(match.arg(equidistant), prob = {
      #cat('Equidistant in probability; ')
      seq.int(from = p.lim[1L], to = p.lim[2L], length.out = p.N)
    }, quantile = {
      #cat('Equidistant in quantile; ')
      if (missing(obs)) stop('Must have `obs`')
      qlim <- quantile(obs, probs = p.lim) # not compute intensive
      ecdf(obs)(seq.int(from = qlim[1L], to = qlim[2L], length.out = p.N))
    }))
  }
  if (!length(p)) stop('will not happen')
  sort.int(p[!duplicated_allequal(p)])
}




