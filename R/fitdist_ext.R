

#' @title Plot \link[fitdistrplus]{fitdist} Object using \CRANpkg{ggplot2}
#' 
#' @description
#' Plot \link[fitdistrplus]{fitdist} object using \CRANpkg{ggplot2}.
#' 
#' @param object a \link[fitdistrplus]{fitdist} object
#' 
#' @param data \link[base]{numeric} or \link[base]{integer} vector, actual observations 
#' 
#' @param type \link[base]{character} scalar, whether the 
#' \code{'density'} (default) or the probability \code{'distribution'} curve should be plotted
#' 
#' @param xlim \link[base]{numeric} vector of length two, the horizontal limit of the figure.
#' Default value is the range of the actual observations.
#' 
#' @param obs.col \link[base]{character} scalar, color to represent the observed values, default black
#' 
#' @param est.col \link[base]{character} scalar, color to represent the estimated values, default red
#' 
#' @param xlab,ylab,title,caption \link[base]{character} scalars, see \link[ggplot2]{labs}
#' 
#' @param ... potential parameters of \link[ggplot2]{stat_function}
#' 
#' @return 
#' 
#' \link{autoplot.fitdist} plots \link[fitdistrplus]{fitdist} object using \CRANpkg{ggplot2}.
#' No value is returned.
#' 
#' @examples
#' 
#' library(fitdistrplus)
#' x1 = rpois(n = 100, lambda = 4)
#' xfit1 = fitdist(x1, distr = 'pois')
#' autoplot(xfit1, xlim = c(-3L, 15L))
#' autoplot(xfit1, type = 'distribution')
#' 
#' x2 = rnorm(n = 1e3L, mean = 2.3, sd = .7)
#' xfit2 = fitdist(x2, distr = 'norm')
#' autoplot(xfit2, type = 'density')
#' autoplot(xfit2, type = 'distribution')
#' 
#' @export
autoplot.fitdist <- function(
  object, data = object[['data']], 
  type = c('density', 'distribution'), 
  xlim = c(min(data), max(data)), 
  obs.col = 'black', est.col = 'red',
  xlab = NULL, ylab = object$distname, title = NULL, caption = NULL,
  ...
) {
  
  if (!length(data)) stop('Provide `data`, or re-run ?fitdistrplus::fitdist with keepdata = TRUE')
  if (anyNA(data)) stop('?fitdistrplus::fitdist will throw error with NA in input `data`')

  if (object$discrete) {
    xlim[1L] <- as.integer(ceiling(xlim[1L]))
    xlim[2L] <- as.integer(floor(xlim[2L]))
    if (xlim[1L] > min(data) || xlim[2L] < max(data)) stop('`xlim` must cover full range of `data`, for discrete distribution')
    xseq <- seq.int(from = xlim[1L], to = xlim[2L], by = 1L) 
  } # else do nothing
  
  type <- match.arg(type)
  
  model_fun <- paste0(switch(type, density = 'd', distribution = 'p'), object$distname)
  
  model_lyr <- if (object$discrete) {
    mp_est <- aes(x = xseq, y = do.call(model_fun, args = c(list(xseq), as.list.default(object$estimate))))
    list(geom_point(mapping = mp_est, colour = est.col), geom_path(mapping = mp_est, colour = est.col, linetype = 2L))
  } else stat_function(fun = model_fun, xlim = xlim, args = as.list.default(object$estimate), colour = est.col, n = 1001L, ...)
  
  data_lyr <- if (object$discrete) {
    y_obs <- tabulate(object$data - xlim[1L] + 1L, nbins = length(xseq)) / length(object$data)
    mp_obs <- aes(x = xseq, y = switch(type, density = y_obs, distribution = cumsum(y_obs))) 
    list(geom_point(mapping = mp_obs, colour = obs.col), geom_path(mapping = mp_obs, colour = obs.col, linetype = 2L))
  } else switch(type, density = geom_histogram(
    mapping = aes_q(x = data, y = quote(..density..)), colour = 'white', alpha = .1, fill = obs.col, bins = 30L
  ), distribution = stat_ecdf(mapping = aes(x = data), geom = 'point', colour = obs.col, shape = 1, size = .5))
  #geom_histogram(mapping = aes_q(x = data, y = quote(cumsum(..density..))), colour = 'white', fill = obs.col, bins = 30L) # why cumulative histgram does not go to 1 ?
  
  ggplot() + data_lyr + model_lyr + 
    labs(x = xlab, y = ylab, title = title, caption = caption) +
    switch(type, distribution = scale_y_continuous(labels = percent)) +
    theme_bw()
  
}


# ?fitdistrplus:::logLik.fitdist not good; wrote to authors
#' @export
logLik.fitdist <- function(object, ...) {
  nobs <- length(data <- object$data)
  if (!nobs) stop('Re-run ?fitdistrplus::fitdist with keepdata = TRUE')
  ret <- object$loglik
  attr(ret, which = 'nobs') <- nobs
  attr(ret, which = 'df') <- length(object$estimate)
  class(ret) <- 'logLik'
  return(ret)
}




