
#' @title Specification of \linkS4class{fmx} Class
#' 
#' @description 
#' Parameters and type of distribution of a one-dimensional finite mixture.
#'
#' @slot distname \link[base]{character} scalar, 
#' name of parametric distribution of the mixture components.
#' Currently, normal (\code{'norm'}) and Tukey's \eqn{g}-&-\eqn{h} (\code{'GH'}) distributions are supported.
#' 
#' @slot pars \link[base]{double} \link[base]{matrix}, 
#' all distribution parameters in the mixture. 
#' Each row corresponds to one component. Each column includes the same parameters of all components.
#' The order of rows corresponds to the (non-strictly) increasing order of the component location parameters.
#' The columns match the formal arguments of the corresponding distribution, 
#' e.g., \code{'mean'} and \code{'sd'} for \link[stats:dnorm]{normal} mixture, 
#' or \code{'A'}, \code{'B'}, \code{'g'} and \code{'h'} for Tukey's \eqn{g}-&-\eqn{h} mixture.
#' 
#' @slot w \link[base]{numeric} \link[base]{vector} of mixing proportions that must sum to 1
#' 
#' @slot data (optional) \link[base]{numeric} \link[base]{vector}, the one-dimensional observations
#' 
#' @slot data.name (optional) \link[base]{character} scalar, a human-friendly name of observations
#'
#' @slot epdf (optional) empirical probability density \link[base]{function} returned by \link[stats]{approxfun}
#' 
# @slot quantile_vv variance-covariance \link[base]{matrix} of selected quantiles 
# (based on the selected probabilities stored in slot \code{@@probs})
#' 
#' @slot vcov_internal (optional) variance-covariance \link[base]{matrix} of the internal (i.e., unconstrained) estimates
#' 
#' @slot vcov (optional) variance-covariance \link[base]{matrix} of the mixture distribution (i.e., constrained) estimates
#' 
# @slot init \linkS4class{fmx} object, the initial values to be sent to \link[stats]{optim}
#' 
#' @slot probs (optional) \link[base]{numeric} vectors of probabilities, where the \link[stats]{quantile}s could be calculated
#' 
# @slot optim a \link[base]{list} returned from \link[stats]{optim}
#' 
#' @export
setClass(Class = 'fmx', slots = c(
  distname = 'character',
  pars = 'matrix',
  w = 'numeric',
  # all below: optional
  data = 'numeric', 
  data.name = 'character',
  epdf = 'function', # prototype of 'function' is `function() NULL`  
  # quantile_vv = 'matrix', # not used
  vcov_internal = 'matrix',
  vcov = 'matrix',
  # init = 'fmx',
  probs = 'numeric'#,
  # optim = 'list' # not used
), prototype = prototype(
  w = 1 # for 1-component
  # data.name = 'observations' # I want to disable this
), validity = function(object) {
  pars <- object@pars
  if (anyNA(pars)) stop('do not allow NA in `fmx` distribution parameter')
  if (is.unsorted(pars[,1L], strictly = FALSE)) { 
    # ?base::is.unsorted.  Note here I use `strictly = FALSE`
    # since log(d) may have large negative value, such that exp(log(d)) may be numerically 0.
    stop('location parameter must be sorted: ', sQuote(pars[,1L])) 
  }
  K <- dim(pars)[1L]
  w <- object@w
  if (K == 1L) {
    if (!isTRUE(all.equal.numeric(w, 1))) stop('slot `w` should be `1` for 1-component') # may not be ?base::identical
  } else {
    if (length(w) != K) stop('slot `w` should be length-K')
    if (!isTRUE(all.equal.numeric(sum(w), 1))) stop('slot `w` must sum up to be `1`') # may not be ?base::identical
  }
  
  # all below: optional
  if (length(object@data)) {
    if (anyNA(object@data)) stop('Observations in \'fmx\' must be free of NA_real_')
    if (object@distname %in% c(distType('continuous'), distType('nonNegContinuous'))) {
      if (!length(body(object@epdf))) stop('empirical pdf must be calculated from the data, for continuous mixture')
    }
  }
}) 



