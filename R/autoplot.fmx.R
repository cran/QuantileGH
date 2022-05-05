

#' @title Plot \linkS4class{fmx} and \linkS4class{fmx_QLMDe} objects using \CRANpkg{ggplot2}
#' 
#' @description 
#' 
#' Plot \linkS4class{fmx} and \linkS4class{fmx_QLMDe} objects using \CRANpkg{ggplot2}.
#' 
#' @param object an \linkS4class{fmx} or \linkS4class{fmx_QLMDe} object
#' 
#' @param type \link[base]{character} scalar.  
#' Option \code{'density'} (default) plots the probability density for \linkS4class{fmx} input or 
#' the histogram and the estimated probability density for \linkS4class{fmx_QLMDe} input.
#' Option \code{'distribution'} plots the cumulative probability distribution for \linkS4class{fmx} input or 
#' the empirical cumulative distribution together with the estimated cumulative distribution function
#' for \linkS4class{fmx_QLMDe} input.
#' 
#' @param data (optional) \link[base]{numeric} vector of the observations.
#' For \linkS4class{fmx_QLMDe} input, the default is \code{object@@data}.
#' 
#' @param epdf ..
#' 
#' @param xlab,ylab,title,caption \link[base]{character} scalars, the 
#' horizontal and vertical label, title and caption.
#' See \code{\link[ggplot2]{xlab}}, \code{\link[ggplot2]{ylab}}, \link[ggplot2]{labs}.
#' 
# @param breaks see parameter \code{breaks} of \link[graphics]{hist.default}, or parameter \code{bins} of \link[ggplot2]{geom_histogram}.
#' 
#' @param hist.col color of the border around the histogram bars, default \code{'white'}.
#' See parameter \code{colour} of \link[ggplot2]{geom_histogram}.
#' 
#' @param hist.fill color of the body of histogram, default \code{'grey95'}.
#' Passed as parameter \code{fill} in \link[ggplot2]{geom_histogram}.
#' 
#' @param curve.col color of the density curve of the fitted finite mixture distribution.
#' Passed as parameter \code{colour} in \link[ggplot2]{stat_function}.
#' 
#' @param xlim horizontal range, see \link[ggplot2]{stat_function}.
#' 
#' @param init \link[base]{logical} scalar, whether to plot the initial estimates used in \code{\link{QLMDe}}, default \code{FALSE}.
#' 
#' @param probs \link[base]{numeric} vector, 
#' the percentages (to be) used in \code{\link{QLMDe}}, can be plotted as vertical lines.
#' Use \code{probs = NULL} to suppress the printing of these lines.
#' 
#' @param origK \link[base]{logical} scalar, 
#' whether to plot the \linkS4class{fmx_QLMDe} at the user-specified number of components \eqn{K}
#' if backward-forward selection on number of component is performed, default \code{FALSE}.
#' 
#' @param ... potential parameters of \link[ggplot2]{stat_function}
#' 
#' @return 
#' 
#' \link{autoplot.fmx} returns a \link[ggplot2]{ggplot} object.
#' 
#' @seealso \link[ggplot2]{autoplot}
#' 
#' @examples 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' curve(dfmx(x, dist = d2), xlim = c(-3, 11))
#' curve(pfmx(x, dist = d2), xlim = c(-3, 11))
#' autoplot(d2)
#' autoplot(d2, type = 'distribution')
#' 
#' @name autoplot_fmx
#' @export
autoplot.fmx <- function(
  object, 
  type = c('density', 'distribution'), 
  data = attr(object, which = 'data', exact = TRUE), # S4 slots, in case those slots does not exists
  epdf = attr(object, which = 'epdf', exact = TRUE),
  probs = attr(object, which = 'p', exact = TRUE),
  init = attr(object, which = 'init', exact = TRUE), 
  origK = attr(object, which = 'orig_K', exact = TRUE),
  xlim = if (!length(data)) qfmx(p = c(.01, .99), dist = object),
  hist.fill = 'grey95', hist.col = 'white',
  curve.col = 1, # black curve, default
  xlab = attr(object, which = 'data.name', exact = TRUE), # S4 slot, 
  ylab = paste(object@distname, 'mixture'), 
  title = TeX(fmx_constraint_brief(object)),
  caption = NULL, #if (nv <- length(v)) paste(nv, 'percentiles to match'),
  ...
) {
  
  # only applicable to continuous distribution
  type <- match.arg(type)
  fun <- switch(type, density = dfmx, distribution = pfmx)
  
  data_lyr <- if (length(data)) {
    switch(type, density = list(
      geom_histogram(mapping = aes_q(x = data, y = quote(..density..)), colour = hist.col, fill = hist.fill, bins = 30L),
      if (is.function(epdf)) stat_function(fun = epdf, n = 1001L, colour = 'grey70', linetype = 2L, ...)
    ), distribution = list(
      stat_ecdf(mapping = aes(x = data), geom = 'step', colour = 'grey70', linetype = 2L)
    ))
  } # else NULL
  
  prob_lyr <- if (length(probs)) {
    v <- if (length(data)) quantile(x = data, probs = probs) else qfmx(p = probs, dist = object)
    vnm <- sprintf('%.1f%%', 1e2*probs)
    list(
      geom_vline(mapping = aes(xintercept = v, colour = vnm), size = .1, show.legend = FALSE),
      geom_text(mapping = aes(x = v, y = 0, colour = vnm, label = vnm, hjust = -.5), size = 3, angle = 90, show.legend = FALSE)
    )
  } # else NULL
  
  ggplot() + data_lyr + prob_lyr + 
    stat_function(fun = fun, args = list(dist = object), xlim = xlim, n = 1001L, colour = curve.col, ...) +
    (if (inherits(init, what = 'fmx')) stat_function(fun = fun, args = list(dist = init), n = 1001L, linetype = 2L, colour = curve.col, ...)) +
    (if (inherits(origK, what = 'fmx')) stat_function(fun = fun, args = list(dist = origK), n = 1001L, colour = curve.col, linetype = 3L, size = .1, ...)) +
    labs(x = xlab, y = ylab, title = title, caption = caption) +
    switch(type, distribution = scale_y_continuous(labels = percent)) +
    theme_bw()
  
}


