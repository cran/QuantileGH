

#' @title Plot \code{'fmx'} and \code{'fmx_QLMDe'} object using vanilla R and/or \pkg{ggplot2}
#' 
#' @description 
#' 
#' Plot \code{'fmx'} and \code{'fmx_QLMDe'} object using vanilla R and/or \pkg{ggplot2}.
#' 
#' @param x,object \code{'fmx'} or \code{'fmx_QLMDe'} object, to match the S3 definition of \code{\link[base]{plot}} and \code{\link[ggplot2]{autoplot}}
#' 
#' @param xlim see \code{\link[graphics]{curve}} and \code{\link[ggplot2]{stat_function}}.
#' 
#' @param xlab,ylab,main,caption,border,breaks see \code{\link[graphics]{hist.default}} and/or \code{\link[ggplot2]{labs}}.
#' 
#' @param hist.col color of the body of histogram, see the parameter \code{col} of \code{\link[graphics]{hist.default}}.
#' 
#' @param curve.col color of the finite mixture density curve, the parameter \code{col} of \code{\link[graphics]{plot.xy}}, or
#' the parameter \code{colour} of \code{\link[ggplot2]{stat_function}}.
#' 
#' @param lty line type of the finite mixture density curve, see \code{\link[graphics]{plot.xy}}, 
#' or parameter \code{linetype} of \code{\link[ggplot2]{stat_function}}.
#' 
#' @param type \code{'character'} value. For an input of \code{'fmx'} object, this argument specifies 
#' whether to generate a plot of probability density function (\code{'pdf'}, default), 
#' or to generate a plot of cumulative distribution function (\code{'cdf'}).
#' For an input of \code{'fmx_QLMDe'} object, this argument specifies 
#' whether to generate a plot of histogram (\code{'hist'}, default) together with the estimated probability density function,
#' or to generate a plot of empirical cumulative distribution function (\code{'ecdf'}) together with the estimated cumulative distribution function.
#' 
#' @param v \code{'double'} vector, the vertical lines to be added to the plot of an \code{'fmx'} object, default \code{double()} indicating no vertical lines.
#' 
#' @param layer \code{'logical'} value, whether to produce the figure as a \pkg{ggplot2} \code{\link[ggplot2]{layer}} (default \code{FALSE}). 
#' 
#' @param n see \code{\link[graphics]{curve}} and/or \code{\link[ggplot2]{stat_function}}.
#' 
#' @param plot_ef,plot_init,plot_p,plot_origK \code{'logical'} value, whether 
#' to plot the empirical probability density function or the empirical probability density function (default \code{TRUE}),
#' to plot the initial estimates (default \code{FALSE}),
#' to plot the percentages at which the sample and theoretical quantiles are matched (default \code{TRUE}),
#' or to plot the \code{'fmx_QLMDe'} at the user-specified number of components \eqn{K} 
#' (if backward-forward selection on number of component is performed; default \code{TRUE}).
#' 
#' @param ... potential parameters of \code{\link[graphics]{curve}} in \code{\link{plot.fmx}}, or 
#' potential parameters of \code{\link[graphics]{hist.default}} in \code{\link{plot.fmx_QLMDe}}.
#' 
#' @return 
#' 
#' \code{\link{plot.fmx}} and \code{\link{plot.fmx_QLMDe}} do not have return value; only \code{\link[graphics]{curve}} is called for plotting.
#' 
#' \code{\link{autoplot.fmx}} and \code{\link{autoplot.fmx_QLMDe}} return \code{'ggplot'} object, created by \code{\link[ggplot2]{ggplot}}.
#' 
#' @name plot_fmx
#' @export
plot.fmx <- function(x, xlim = qfmx(p = c(.01, .99), dist = dist), ylab = paste('Finite Mixture of', dist@distname), n = 501L, ...) {
  dist <- x
  curve(dfmx(x, dist = dist), xlim = xlim, ylab = ylab, n = n, ...) # `x` is simply the 1st arg in ?graphics::curve
}


# ?ggplot2::autoplot
#' @rdname plot_fmx
#' @export
autoplot.fmx <- function(
  object, type = c('pdf', 'cdf'), xlim = qfmx(p = c(.005, .995), dist = object),
  v = double(), # default is no vertical lines at all
  layer = FALSE,
  lty = 1L, curve.col = 1L,
  ylab = paste(object@distname, 'mixture'), 
  n = 501L,
  ...
) {
  
  #if (any(object@distname == .distr_p_int)) stop('autoplot only applicable to continuous distribution')
  type <- match.arg(type)
  
  if (length(v)) {
    if (!is.double(v) || anyNA(v)) stop('do not allow NA in `v`')
    v <- v[v > xlim[1L] & v < xlim[2L]] # vertical lines outside specified `xlim` are removed
  }
  if (has_v <- length(v)) {
    vnm <- names(v)
    if (!length(vnm)) vnm <- sprintf('%.1f', v)
  }
  
  fun <- switch(type, pdf = dfmx, cdf = pfmx, stop(sQuote(type), ' not programed yet'))
  lyr <- list(
    stat_function(fun = fun, xlim = xlim, args = list(dist = object), n = n, linetype = lty, colour = curve.col),
    (if (has_v) geom_vline(mapping = aes(xintercept = v, colour = vnm), size = .4, alpha = .3, show.legend = FALSE)),
    (if (has_v) geom_text(mapping = aes(x = v, y = 0, colour = vnm, label = vnm, hjust = -.5), angle = 90, show.legend = FALSE))
  )
  if (layer) return(lyr)
  
  ggplot() + lyr +
    labs(y = ylab, caption = if (nv <- length(v)) paste(nv, 'percentiles to match')) +
    (if (type == 'cdf') scale_y_continuous(labels = percent)) +
    theme(panel.background = element_rect(fill = 'transparent', colour = 'black'),
          panel.grid.major = element_line(colour = 'grey95', size = .2), panel.grid.minor = element_blank(),
          axis.title.x = element_blank())
}







#' @rdname plot_fmx
#' @export
plot.fmx_QLMDe <- function(
  x, 
  xlab = x@data.name, main = NULL, 
  hist.col = 'grey90', breaks = 'Sturges', border = 'white',
  curve.col = 'black', lty = 1L, n = 501L,
  ...
) {
  hist.default(x@data, freq = FALSE, col = hist.col, breaks = breaks, border = border, xlab = xlab, main = main, ...)
  fitted <- x
  # `x` is simply the 1st arg in ?graphics::curve
  curve(fitted@epdf(x), lty = 2L, add = TRUE, n = n)
  curve(dfmx(x, dist = fitted), col = curve.col, lty = lty, add = TRUE, n = n, ...)
}



# ?ggplot2::autoplot
#' @rdname plot_fmx
#' @export
autoplot.fmx_QLMDe <- function(
  object, type = c('hist', 'ecdf'),
  plot_ef = TRUE, plot_init = FALSE, plot_p = TRUE, plot_origK = TRUE,
  curve.col = 1, # black curve, default
  xlim = c(min(object@data), max(object@data)),
  xlab = object@data.name, 
  main = NULL,
  caption = if (plot_p) paste(length(object@p), 'percentiles matched'),
  layer = FALSE,
  n = 501L,
  ...
) {
  
  #if (any(object@distname == .distr_p_int)) stop('autoplot only applicable to continuous distribution')
  type <- match.arg(type)
  type1 <- switch(type, hist = 'pdf', ecdf = 'cdf')
  
  lyr <- list(
    if (plot_ef) stat_function(fun = switch(type, hist = object@epdf, ecdf = ecdf(object@data)), xlim = xlim, n = n, colour = 'grey70', linetype = 2L, size = 1.5),
    if (plot_init) autoplot.fmx(object@init, type = type1, xlim = xlim, n = n, colour = curve.col, linetype = 2L, layer = TRUE, ...), 
    autoplot.fmx(object, type = type1, v = if (plot_p) quantile(x = object@data, probs = object@p), xlim = xlim, n = n,  colour = curve.col, layer = TRUE, ...),
    if (plot_origK && length(y_origK <- attr(object, which = 'orig_K', exact = TRUE))) autoplot.fmx(y_origK, type = type1, xlim = xlim, n = n, colour = curve.col, linetype = 3L, size = 1.5, layer = TRUE, ...)
  )
  if (layer) return(lyr)
  
  ggplot() + 
    switch(type, hist = geom_histogram(mapping = aes_q(x = object@data, y = quote(..density..)), colour = 'white', alpha = .1, bins = 30L)) +
    lyr + 
    labs(x = xlab, y = ylab, main = main, caption = caption) +
    (if (type == 'ecdf') scale_y_continuous(labels = percent)) +
    theme(panel.background = element_rect(fill = 'transparent', colour = 'black'),
          panel.grid.major = element_line(colour = 'grey95', size = .2), panel.grid.minor = element_blank(),
          axis.title.y = element_blank())
}

