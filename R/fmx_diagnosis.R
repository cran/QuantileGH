
#' @title Diagnoses for \linkS4class{fmx} Estimates
#' 
#' @description 
#' 
#' Diagnoses for \linkS4class{fmx} estimates.
#' 
#' @param object \linkS4class{fmx} object, or an R object convertible to an \linkS4class{fmx} object
#' 
#' @param data \link[base]{double} \link[base]{vector}, observed data.
#' Default is `object@@data`, the data used for estimation.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' Function [Kolmogorov_fmx()] calculates Kolmogorov distance.
#' 
#' @returns
#' Functions [Kolmogorov_fmx()], [KullbackLeibler_fmx()], [CramerVonMises_fmx()] 
#' all return \link[base]{numeric} scalars.
#' 
#' @name fmx_diagnosis
#' @export
Kolmogorov_fmx <- function(object, data = object@data, ...) {
  
  if (!length(object)) return(NA_real_)
  if (length(ret <- attr(object, which = 'Kolmogorov'))) return(ret) # do not use `@` for back compatibility
  
  object <- as.fmx(object) # uses S3
  Kolmogorov_dist(data, null = pfmx, dist = object)
  
}






#' @details
#' Function [KullbackLeibler_fmx()] calculates Kullback-Leibler divergence.
#' The R code is adapted from `LaplacesDemon::KLD`.
#' 
#' @rdname fmx_diagnosis
#' @export
KullbackLeibler_fmx <- function(object, data = object@data, ...) {
  
  if (!length(object)) return(NA_real_)
  if (length(ret <- attr(object, which = 'KullbackLeibler'))) return(ret) # do not use `@` for back compatibility
  
  object <- as.fmx(object) # uses S3
  if (!length(data)) stop('must provide actual observations')
  
  # read ?LaplacesDemon::KLD carefully
  # `px` is model-based; `py` is empirical
  px <- dfmx(data, dist = object, log = FALSE)
  py <- if (object@distname %in% c(distType('continuous'), distType('nonNegContinuous'))) {
    object@epdf(data)
  } else (tabulate(data, nbins = max(data)) / length(data))[data]
  if (any(!is.finite(px), !is.finite(py))) stop('`px` and `py` must have finite values.')
  px <- pmax.int(px, .Machine$double.xmin)
  py <- pmax.int(py, .Machine$double.xmin)
  px <- px/sum(px) # um..
  py <- py/sum(py) # normalization ..
  return(sum(py * (log(py) - log(px))))
  
}







#' @details
#' Function [CramerVonMises_fmx()] calculates Cramer-von Mises quadratic distance 
#' (via \link[goftest]{cvm.test}).
#' 
#' @seealso `dgof::cvmf.test`
#' @importFrom goftest cvm.test
#' @rdname fmx_diagnosis
#' @export
CramerVonMises_fmx <- function(object, data = object@data, ...) {
  
  if (!length(object)) return(NA_real_)
  if (length(ret <- attr(object, which = 'CramerVonMises'))) return(ret) # do not use `@` for back compatibility
  
  object <- as.fmx(object) # uses S3
  if (!length(data)) stop('must provide actual observations')
  
  # cvm.test(x = data, null = pfmx, dist = object, nullname = nullname, ...) # cannot deal complex parameters such as my 'fmx'
  unname(cvm.test(x = data, null = function(q) {
    pmin.int(pmax.int(pfmx(q, dist = object), .Machine$double.eps), 1 - .Machine$double.eps)
    #pfmx(q, dist = object)
  }, estimated = TRUE, nullname = '', ...)[['statistic']])
  # why each time different???

}




