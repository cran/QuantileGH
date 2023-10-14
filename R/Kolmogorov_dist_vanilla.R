


#' @title One-Sample Kolmogorov Distance
#' 
#' @description 
#' 
#' To calculate the one-sample Kolmogorov distance between observations and 
#' a distribution. 
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, observations \eqn{x}
#' 
#' @param null cumulative distribution \link[base]{function}
#' 
#' @param alternative \link[base]{character} scalar,
#' alternative hypothesis, either `'two.sided'` (default), `'less'`, or `'greater'`
#' 
#' @param ... additional arguments of `null`
#' 
#' @returns 
#' Function [Kolmogorov_dist()] returns a \link[base]{numeric} scalar.
#' 
#' @details 
#' Function [Kolmogorov_dist()] is different from \link[stats]{ks.test} in the
#' following aspects
#' \itemize{
#' \item {Ties in observations are supported.  
#' The step function of empirical distribution is inspired by \link[stats]{ecdf}.
#' This is superior than `(0:(n - 1))/n` in \link[stats]{ks.test}.}
#' \item {Discrete distribution (with discrete observation) is supported.}
#' \item {Discrete distribution (with continuous observation) is not supported yet.
#' This will be an easy modification in future.}
#' \item {Only the one-sample Kolmogorov distance, not the one-sample Kolmogorov test, is returned, 
#' to speed up the calculation.}
#' }
#' 
#' @seealso 
#' \link[stats]{ks.test}
#' 
#' @examples 
#' # from ?stats::ks.test
#' x1 = rnorm(50)
#' ks.test(x1+2, y = pgamma, shape = 3, rate = 2)
#' Kolmogorov_dist(x1+2, null = pgamma, shape = 3, rate = 2) # exactly the same
#' 
#' # discrete distribution
#' x2 <- rnbinom(n = 1e2L, size = 500, prob = .4)
#' suppressWarnings(ks.test(x2, y = pnbinom, size = 500, prob = .4)) # warning on ties
#' Kolmogorov_dist(x2, null = pnbinom, size = 500, prob = .4) # wont be the same
#' 
#' @export
Kolmogorov_dist <- function(x, null, alternative = c('two.sided', 'less', 'greater'), ...) {
  
  if (!length(x) || !is.numeric(x) || anyNA(x)) stop('input data `x` must be numeric, free of missingness')
  if (is.character(null)) 
    null <- get(null, mode = 'function', envir = parent.frame())
  if (!is.function(null)) stop('`null` must be a function')
  
  n <- length(x)
  if (n < 1L) stop('not enough observations')
  
  x <- sort.int(x)
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
