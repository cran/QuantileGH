


#' @title One-Sample Kolmogorov Distance
#' 
#' @description 
#' 
#' To calculate the one-sample Kolmogorov distance between observations and 
#' a distribution. 
#' 
#' @param obs \link[base]{numeric} or \link[base]{integer} \link[base]{vector}, observations
#' 
#' @param null cumulative distribution \link[base]{function}
#' 
#' @param alternative \link[base]{character} scalar,
#' alternative hypothesis, either \code{'two.sided'} (default), \code{'less'}, or \code{'greater'}
#' 
#' @param ... additional arguments of \code{null}
#' 
#' @return 
#' \link{Kolmogorov_dist} returns a \link[base]{numeric} scalar.
#' 
#' @details 
#' \link{Kolmogorov_dist} is different from \link[stats]{ks.test} in the
#' following aspects
#' \itemize{
#' \item {Ties in observations are supported.  
#' The step function of empirical distribution is inspired by \link[stats]{ecdf}.
#' This is superior than \code{(0:(n - 1))/n} of \code{stats:::ks.test.default}.}
#' \item {Discrete distribution (with discrete observation) is supported.}
#' \item {Discrete distribution (with continuous observation) is not supported yet.
#' This will be an easy modification in future.}
#' \item {Only the one-sample Kolmogorov distance, not the one-sample Kolmogorov test, is returned, 
#' to speed up the calculation.}
#' }
#' 
#' @seealso \link[stats]{ks.test}
#' 
#' @examples 
#' # from ?stats::ks.test
#' x1 <- rnorm(50)
#' ks.test(x1+2, y = pgamma, shape = 3, rate = 2)
#' Kolmogorov_dist(x1+2, null = pgamma, shape = 3, rate = 2) # exactly the same
#' 
#' # discrete distribution
#' x2 <- rnbinom(n = 1e2L, size = 500, prob = .4)
#' suppressWarnings(ks.test(x2, y = pnbinom, size = 500, prob = .4)) # warning on ties
#' Kolmogorov_dist(x2, null = pnbinom, size = 500, prob = .4) # wont be the same
#' 
#' @export
Kolmogorov_dist <- function(obs, null, alternative = c('two.sided', 'less', 'greater'), ...) {
  
  if (!is.numeric(obs) || anyNA(obs)) stop('input must be numeric, free of missingness')
  if (is.character(null)) 
    null <- get(null, mode = 'function', envir = parent.frame())
  if (!is.function(null)) stop('`null` must be a function')
  
  n <- length(obs)
  if (n < 1L) stop('not enough `obs` data')
  
  x <- sort.int(obs)
  ux <- unique.default(x)
  epdf <- cumsum(tabulate(match(x, table = ux)))/n - 1/n # `-1/n` to match the behavior of ?stats::ks.test
  
  ret <- null(ux, ...) - epdf
  # only `$statistic` from ?stats:::ks.test.default
  switch(match.arg(alternative), two.sided = {
    max(c(ret, 1/n - ret))
  }, greater = {
    max(1/n - ret)
  }, less = {
    max(ret)
  })
  
}
