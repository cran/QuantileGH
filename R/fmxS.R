
#' @title \linkS4class{fmxS}: Multiple \linkS4class{fmx} objects
#' 
#' @description
#' ..
#' 
#' @slot .Data \link[base]{list} of \linkS4class{fmx} objects
#' 
#' @slot data \link[base]{numeric} \link[base]{vector}
#' 
#' @slot data.name \link[base]{character} scalar
#' 
#' @name fmxS
#' @aliases fmxS-class
#' @export
setClass(Class = 'fmxS', contains = 'list', slots = c(
  data = 'numeric',
  data.name = 'character'
))


#' @rdname fmxS
#' 
#' @param ... multiple \linkS4class{fmx} objects, or objects convertible
#' to \linkS4class{fmx} class via function [as.fmx()]
#' 
#' @returns
#' Function [fmxS()] returns an \linkS4class{fmxS} object.
#' 
#' @examples 
#' 
#' library(fitdistrplus)
#' set.seed(1234); x = rnorm(n = 1e3L)
#' f1 = fitdist(x, distr = 'norm')
#' f2 = fitdist(x, distr = 'GH', start = as.list.default(letterValue(x)))
#' aa = fmxS(a = f1, b = f2)
#' summary(aa)
#' autoplot(aa, type = 'density')
#' autoplot(aa, type = 'distribution')
#' 
#' a1 = fmx('GH', A = c(7,9), B = c(.8, 1.2), g = c(.3, 0), h = c(0, .1), w = c(1, 1))
#' a2 = fmx('GH', A = c(6,9), B = c(.8, 1.2), g = c(-.3, 0), h = c(.2, .1), w = c(4, 6))
#' a = fmxS(a1, a2)
#' (p = autoplot(a, type = 'distribution') + coord_flip())
#' p + labs(x = 'new xlab', y = 'new ylab')
#' p + theme(legend.position = 'none')
#' 
#' @export
fmxS <- function(...) {
  dots <- list(...)
  dots <- dots[!duplicated.default(dots)] # ?base::unique.default drops names
  dots <- lapply(dots, FUN = as.fmx)

  if (is.null(nm <- names(dots))) names(dots) <- gsub('^GH2\\: ', replacement = '', vapply(dots, FUN = getTeX, FUN.VALUE = ''))
  if (anyNA(nm) || !all(nzchar(nm)) || anyDuplicated.default(nm)) stop('names of `dots` must be unique')
  
  if ((length(datas <- unique.default(lapply(dots, FUN = slot, name = 'data'))) != 1L) ||
      (length(datanms <- unique.default(lapply(dots, FUN = slot, name = 'data.name'))) != 1L)) stop('`dots` are not based on the same data')
  data <- datas[[1L]]
  data.name <- datanms[[1L]]
  
  new(Class = 'fmxS', dots, 
      data = data, data.name = data.name)
}





#' @title Plot \linkS4class{fmxS} Object
#' 
#' @description
#' ..
#' 
#' @param object \linkS4class{fmxS} object
#' 
#' @param type \link[base]{character} scalar, `'density'` (default) or `'distribution'`
#' 
#' @param xlim ..
#' 
#' @param ... ..
#' 
#' @returns
#' Function [autoplot.fmxS()] returns a \link[ggplot2]{ggplot} figure.
#' 
#' @export
autoplot.fmxS <- function(
    object, 
    type = c('density', 'distribution'), 
    xlim = if (length(data)) range.default(data) else range.default(lapply(dots, FUN = qfmx, p = c(.01, .99))),
    ...
) {
  
  N <- length(data <- object@data)
  data.name <- object@data.name
  dots <- asS3(object)
  
  type <- match.arg(type)
  
  args <- lapply(dots, FUN = function(i) list(dist = i))
  
  if (N) {

    epdf_continuous <- list(
      geom_histogram(mapping = aes(x = data, y = after_stat(density)), colour = 'white', fill = 'grey90', bins = 30L), # creates ylab 'density'
      geom_density(mapping = aes(x = data), colour = 'grey70', lty = 4L, size = 1.2)  # creates ylab 'y'
    )
    ecdf_continuous <- stat_ecdf(mapping = aes(x = data), geom = 'step', pad = FALSE)
    pdf_continuous <- geom_function_args(fun = dfmx, args = args, ...) # , alpha = .7 
    cdf_continuous <- geom_function_args(fun = pfmx, args = args, ...) # , alpha = .7 
    
    # fig_message <- .cyan(c('Figure of ', sQuote(data.name), ' created'))
    xlab <- sprintf(fmt = '%s (n=%d)', data.name, N)
    
  } else { # no data
    
    epdf_continuous <- ecdf_continuous <- xlab <- NULL
    # xlim(xlim), 
    # ?ggplot2::xlim must, default `xlim` determined by ?ggplot2::geom_histogram and ?ggplot2::stat_ecdf can be different
    # masked for now, as I do not know the true `xlim` determined by ?ggplot2::geom_histogram
    # let user set xlim.
    pdf_continuous <- geom_function_args(fun = dfmx, args = args, xlim = xlim, ...)
    cdf_continuous <- geom_function_args(fun = pfmx, args = args, xlim = xlim, ...)
    
  }
  
  p <- ggplot() + switch(type, density = list(
    epdf_continuous, pdf_continuous
  ), distribution = list(
    ecdf_continuous, cdf_continuous
  )) + labs(x = xlab, y = NULL, colour = 'Models')
  
  #if (length(msg <- fig_message)) on.exit(message(msg))
  return(p)
}


