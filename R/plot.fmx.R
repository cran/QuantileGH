

#' @title Plot \code{\linkS4class{fmx}} and \code{\linkS4class{fmx_QLMDe}} object using vanilla R and/or \pkg{ggplot2}
#' 
#' @description 
#' 
#' Plot \code{\linkS4class{fmx}} and \code{\linkS4class{fmx_QLMDe}} object using vanilla R and/or \pkg{ggplot2}.
#' 
#' @param x,object \code{\linkS4class{fmx}} or \code{\linkS4class{fmx_QLMDe}} object, to match the S3 definition of \code{\link[base]{plot}} and \code{\link[ggplot2]{autoplot}}
#' 
#' @param xlab,ylab,main,title,caption see \code{\link[graphics]{curve}}, \code{\link[graphics]{hist.default}} and/or \code{\link[ggplot2]{labs}}.
#' 
#' @param breaks see parameter \code{breaks} of \code{\link[graphics]{hist.default}}, or parameter \code{bins} of \code{\link[ggplot2]{geom_histogram}}.
#' 
#' @param border color of the border around the histogram bars, default \code{'white'}.
#' See parameter \code{border} of \code{\link[graphics]{hist.default}}, 
#' or parameter \code{colour} of \code{\link[ggplot2]{geom_histogram}}.
#' 
#' @param hist.col color of the body of histogram, default \code{'grey95'}.
#' See parameter \code{col} of \code{\link[graphics]{hist.default}}, 
#' or parameter \code{fill} of \code{\link[ggplot2]{geom_histogram}}.
#' 
#' @param curve.col color of the density curve of the fitted finite mixture distribution.
#' See parameter \code{col} of \code{\link[graphics]{plot.xy}},
#' or parameter \code{colour} of \code{\link[ggplot2]{stat_function}}.
#' 
#' @param lty line type of the finite mixture density curve, 
#' see parameter \code{lty} of \code{\link[graphics]{plot.xy}}, 
#' or parameter \code{linetype} of \code{\link[ggplot2]{stat_function}}.
#' 
#' @param type \code{\link[base]{character}} value. For an input of \code{\linkS4class{fmx}} object, this argument specifies 
#' whether to generate a plot of probability density function (\code{'pdf'}, default), 
#' log of probability density function (\code{'log-pdf'}),
#' or cumulative distribution function (\code{'cdf'}).
#' For an input of \code{\linkS4class{fmx_QLMDe}} object, this argument specifies 
#' whether to generate a plot of histogram (\code{'hist'}, default) together with the estimated probability density function,
#' or to generate a plot of empirical cumulative distribution function (\code{'ecdf'}) together with the estimated cumulative distribution function.
#' 
#' @param v \code{\link[base]{double}} vector, the vertical lines to be added to the plot of an \code{\linkS4class{fmx}} object, default \code{double()} indicating no vertical lines.
#' 
#' @param layer \code{\link[base]{logical}} value, whether to produce the figure as a \code{\link[ggplot2]{layer}} (default \code{FALSE}). 
#' 
#' @param xlim,n see \code{\link[graphics]{curve}} and/or \code{\link[ggplot2]{stat_function}}.
#' 
#' @param plot_ef \code{\link[base]{logical}} value, whether to plot 
#' the empirical probability density or cumulative probability function, default \code{TRUE}.
#' 
#' @param plot_init \code{\link[base]{logical}} value, whether to plot the initial estimates used in \code{\link{QLMDe}}, default \code{FALSE}.
#' 
#' @param plot_p \code{\link[base]{logical}} value, whether to plot the percentages used in \code{\link{QLMDe}}, default \code{TRUE}.
#' 
#' @param plot_origK \code{\link[base]{logical}} value, 
#' whether to plot the \code{\linkS4class{fmx_QLMDe}} at the user-specified number of components \eqn{K}
#' if backward-forward selection on number of component is performed, default \code{FALSE}.
#' 
#' @param ... potential parameters of \code{\link[graphics]{curve}} in \code{\link{plot.fmx}}, 
#' potential parameters of \code{\link[graphics]{hist.default}} in \code{\link{plot.fmx_QLMDe}}, or
#' potential parameters of \code{\link[ggplot2]{stat_function}} in \code{\link{autoplot.fmx}} and \code{\link{autoplot.fmx_QLMDe}}.
#' 
#' @return 
#' 
#' \code{\link{plot.fmx}} and \code{\link{plot.fmx_QLMDe}} do not have return value; only \code{\link[graphics]{curve}} is called for plotting.
#' 
#' \code{\link{autoplot.fmx}} and \code{\link{autoplot.fmx_QLMDe}} return \code{'ggplot'} object, created by \code{\link[ggplot2]{ggplot}}.
#' 
#' @seealso \code{\link[graphics]{plot}}, \code{\link[ggplot2]{autoplot}}
#' 
#' @examples 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' plot(d2)
#' plot(d2, type = 'cdf')
#' autoplot(d2)
#' autoplot(d2, type = 'cdf')
#' 
#' @name plot_fmx
#' @export
plot.fmx <- function(
  x, type = c('pdf', 'cdf', 'log_pdf'),
  xlim = qfmx(p = c(.01, .99), dist = dist), n = 501L, 
  xlab = '', # must use `''` to suppress the label, \code{NULL} wont work
  ylab = paste(dist@distname, 'mixture'), 
  ...
) {
  dist <- x
  # cannot use ?base::match.arg inside ?graphics::curve
  switch(match.arg(type), pdf = {
    curve(dfmx(x, dist = dist), xlim = xlim, xlab = xlab, ylab = ylab, n = n, ...)
  }, log_pdf = {
    curve(dfmx(x, dist = dist, log = TRUE), xlim = xlim, xlab = xlab, ylab = ylab, n = n, ...)
  }, cdf = {
    curve(pfmx(x, dist = dist), xlim = xlim, xlab = xlab, ylab = ylab, n = n, ...) 
  })
}



#' @rdname plot_fmx
#' @export
autoplot.fmx <- function(
  object, type = c('pdf', 'cdf'), xlim = qfmx(p = c(.01, .99), dist = object),
  v = double(), # default is no vertical lines at all
  layer = FALSE,
  lty = 1L, curve.col = 1L,
  xlab = NULL, ylab = paste(object@distname, 'mixture'), title = NULL,
  caption = NULL, #if (nv <- length(v)) paste(nv, 'percentiles to match'),
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
    stat_function(fun = fun, xlim = xlim, args = list(dist = object), n = n, linetype = lty, colour = curve.col, ...),
    (if (has_v) geom_vline(mapping = aes(xintercept = v, colour = vnm), size = .1, show.legend = FALSE)), # , alpha = .3
    (if (has_v) geom_text(mapping = aes(x = v, y = 0, colour = vnm, label = vnm, hjust = -.5), size = 3, angle = 90, show.legend = FALSE))
  )
  if (layer) return(lyr)
  
  ggplot() + lyr +
    labs(x = xlab, y = ylab, title = title, caption = caption) +
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
  hist.col = 'grey95', breaks = 'Sturges', border = 'white',
  curve.col = 'black', lty = 1L, n = 501L,
  ...
) {
  hist.default(x@data, freq = FALSE, col = hist.col, breaks = breaks, border = border, xlab = xlab, main = main, ...)
  fitted <- x
  # `x` is simply the 1st arg in ?graphics::curve
  curve(fitted@epdf(x), col = 'grey70', lty = 2L, add = TRUE, n = n)
  curve(dfmx(x, dist = fitted), col = curve.col, lty = lty, add = TRUE, n = n, ...)
}




#' @rdname plot_fmx
#' @export
autoplot.fmx_QLMDe <- function(
  object, type = c('hist', 'ecdf'),
  plot_ef = TRUE, plot_init = FALSE, plot_p = TRUE, plot_origK = FALSE,
  hist.col = 'grey95', border = 'white',
  curve.col = 1, # black curve, default
  xlim = c(min(object@data), max(object@data)),
  xlab = object@data.name, ylab = paste(object@distname, 'mixture'), title = NULL,
  caption = if (plot_p) paste(length(object@p), 'percentiles matched'),
  layer = FALSE,
  n = 501L,
  ...
) {
  
  #if (any(object@distname == .distr_p_int)) stop('autoplot only applicable to continuous distribution')
  type <- match.arg(type)
  type1 <- switch(type, hist = 'pdf', ecdf = 'cdf')
  
  lyr <- list(
    if (plot_ef) stat_function(fun = switch(type, hist = object@epdf, ecdf = ecdf(object@data)), xlim = xlim, n = n, colour = 'grey70', linetype = 2L, ...),
    if (plot_init) autoplot.fmx(object@init, type = type1, xlim = xlim, n = n, colour = curve.col, linetype = 2L, layer = TRUE, ...), 
    autoplot.fmx(object, type = type1, v = if (plot_p) setNames(quantile(x = object@data, probs = object@p), nm = sprintf('%.1f%%', 1e2*object@p)), 
                 xlim = xlim, n = n, curve.col = curve.col, layer = TRUE, ...),
    if (plot_origK && length(y_origK <- attr(object, which = 'orig_K', exact = TRUE))) autoplot.fmx(y_origK, type = type1, xlim = xlim, n = n, colour = curve.col, lty = 3L, size = 1.5, layer = TRUE, ...)
  )
  if (layer) return(lyr)
  
  data_lyr <- switch(type, hist = {
    geom_histogram(mapping = aes_q(x = object@data, y = quote(..density..)), colour = border, fill = hist.col, bins = 30L)
  }, ecdf = {
    stat_ecdf(mapping = aes(x = object@data), geom = 'step', colour = hist.col)
  })
  
  ggplot() + data_lyr + lyr + 
    labs(x = xlab, y = ylab, title = title, caption = caption) +
    (if (type == 'ecdf') scale_y_continuous(labels = percent)) +
    theme(panel.background = element_rect(fill = 'transparent', colour = 'black'),
          panel.grid.major = element_line(colour = 'grey95', size = .2), panel.grid.minor = element_blank())
}



#' @export
autoplot.fitdist <- function(
  object, type = c('hist', 'ecdf'), 
  xlim = c(min(data), max(data)), 
  hist.col = 'grey95', border = 'white',
  curve.col = 1, # black curve, default
  xlab = NULL, ylab = object$distname, title = NULL, caption = NULL,
  layer = FALSE, n = 501L, 
  ...
) {
  if (!length(data <- object$data)) stop('Re-run ?fitdistrplus::fitdist with keepdata = TRUE')
  if (anyNA(data)) stop('?fitdistrplus::fitdist will throw errow with NA in input `data`')
  if (object$discrete) stop('only continuous \'fitdist\' supported, yet')
  type <- match.arg(type)
  lyr <- list(
    stat_function(
      fun = paste0(switch(type, hist = 'd', ecdf = 'p'), object$distname), 
      xlim = xlim, args = as.list.default(object$estimate), colour = curve.col, n = n, ...)
  )
  if (layer) return(lyr)
  
  data_lyr <- switch(type, hist = {
    geom_histogram(mapping = aes_q(x = data, y = quote(..density..)), colour = border, fill = hist.col, bins = 30L)
  }, ecdf = {
    # why cumulative histgram does not go to 1 ?
    #geom_histogram(mapping = aes_q(x = data, y = quote(cumsum(..density..))), colour = border, fill = hist.col, bins = 30L)
  })
  
  ggplot() + data_lyr + lyr + 
    labs(x = xlab, y = ylab, title = title, caption = caption) +
    (if (type == 'ecdf') scale_y_continuous(labels = percent)) +
    theme(panel.background = element_rect(fill = 'transparent', colour = 'black'),
          panel.grid.major = element_line(colour = 'grey95', size = .2), panel.grid.minor = element_blank())
}




