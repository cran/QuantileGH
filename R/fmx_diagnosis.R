



#' @title Diagnoses for \linkS4class{fmx} Estimates
#' 
#' @description 
#' 
#' Diagnoses for \linkS4class{fmx} estimates.
#' 
#' @param x \linkS4class{fmx} object, or an R object convertible to an \linkS4class{fmx} object
#' 
#' @param data \link[base]{double} \link[base]{vector}, the actual observations, default value is \code{@@data}
#' 
#' @param type \link[base]{character} scalar, currently supports 
#' \code{'Kolmogorov'} distance (default),
#' \code{'KullbackLeibler'} divergence,
#' and \code{'CramerVonMises'} test of goodness-of-fit
#' 
#' @param nullname \link[base]{character} scalar, see \link[goftest]{cvm.test}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @return 
#' 
#' \link{fmx_diagnosis} returns either
#' 
#' \itemize{
#' 
#' \item {a \link[base]{numeric} scalar of \code{'Kolmogorov'} distance;}
#' 
#' \item {a \link[base]{list} when \code{type = 'KullbackLeibler'}, which is returned from \code{LaplacesDemon::KLD}.}
#' 
#' \item {an \link[goftest:cvm.test]{htest} object when \code{type = 'CramerVonMises'},
#'  in which the element \code{$statistic} is the Cramer-Von Mises quadratic distance.}
#' 
#' }
#' 
#' @seealso 
#' \code{LaplacesDemon::KLD} \link[goftest]{cvm.test} 
#' \link[stats]{ks.test}
#' 
#' @note \link[dgof]{cvm.test}
#' 
#' @importFrom goftest cvm.test
#' 
#' @export
fmx_diagnosis <- function(
    x, data = x@data,
    type = c('Kolmogorov', 'KullbackLeibler', 'CramerVonMises'), 
    nullname = deparse1(substitute(x)), 
    ...
) {
  
  force(nullname)
  
  x <- as.fmx(x) # uses S3
  if (!length(data)) stop('must provide actual observations')
  
  switch(match.arg(type), Kolmogorov = {
    
    #suppressWarnings(ks.test(x = data, y = pfmx, dist = x, ...))# Warnings on ties
    return(Kolmogorov_dist(obs = data, null = pfmx, dist = x))
    
  }, KullbackLeibler = {
    
    # read ?LaplacesDemon::KLD carefully
    # `px` is model-based; `py` is empirical
    px <- dfmx(data, dist = x, log = FALSE)
    py <- if (x@distname %in% c(distType('continuous'), distType('nonNegContinuous'))) {
      x@epdf(data)
    } else (tabulate(data, nbins = max(data)) / length(data))[data]
    if (any(!is.finite(px), !is.finite(py))) stop('`px` and `py` must have finite values.')
    px <- pmax.int(px, .Machine$double.xmin)
    py <- pmax.int(py, .Machine$double.xmin)
    px <- px/sum(px) # um..
    py <- py/sum(py) # normalization ..
    return(sum(py * (log(py) - log(px))))
    
  }, CramerVonMises = {
    
    # cvm.test(x = data, null = pfmx, dist = x, nullname = nullname, ...) # cannot deal complex parameters such as my 'fmx'
    cvm.test(x = data, null = function(q) {
      pmin.int(pmax.int(pfmx(q, dist = x), .Machine$double.eps), 1 - .Machine$double.eps)
      #pfmx(q, dist = x)
    }, estimated = TRUE, nullname = nullname, ...)
    
  })
  
}




