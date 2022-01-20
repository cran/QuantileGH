
#' @title Specification of \code{\linkS4class{fmx}} Class
#' 
#' @description 
#' Parameter specification for a one-dimensional finite mixture distribution.
#'
#' @slot distname \code{\link[base]{character}} value for the name of parametric distribution of mixture components.
#' Currently normal (\code{'norm'}) and Tukey's \eqn{g}-&-\eqn{h} (\code{'GH'}) distributions are supported.
#' 
#' @slot parM \code{\link[base]{matrix}} of all distribution parameters in the mixture. 
#' Each row corresponds to one component. Each column includes the same parameters of all components.
#' The order of rows corresponds to the (non-strictly) increasing order of the component location parameters.
#' The column names, e.g., \code{mean} and \code{sd} for norm distribution (see \code{\link[stats]{dnorm}}), 
#' or \code{A}, \code{B}, \code{g} and \code{h} for Tukey's \eqn{g}-&-\eqn{h} distribution (see \code{\link{dGH}}).
#' 
#' @slot w \code{\link[base]{numeric}} vector of mixing proportions that must sum to 1
#' 
# @details ..
#' 
#' @examples 
#' ?`fmx-class`
#' 
#' @export
setClass(Class = 'fmx', slots = c(
  distname = 'character',
  parM = 'matrix',
  w = 'numeric'
), prototype = prototype(
  w = 1 # for 1-component
), validity = function(object) {
  parM <- object@parM
  if (anyNA(parM)) stop('do not allow NA in `fmx` distribution parameter')
  if (.Internal(is.unsorted(parM[,1L], FALSE))) stop('location parameter must be sorted') 
  K <- dim(parM)[1L]
  w <- object@w
  if (K == 1L) {
    #if (!identical(w, 1)) stop('slot `w` should be `1` for 1-component')
    if (!isTRUE(all.equal.numeric(w, 1))) stop('slot `w` should be `1` for 1-component') # may not be ?base::identical
  } else if (length(w) != K) stop('slot `w` should be length-K')
  # ?base::is.unsorted.  Note here I use `strictly = FALSE`, 
  # since log(d) may have large negative value, such that exp(log(d)) may be numerically 0.
}) 



#' @title Specification of \code{\linkS4class{fmx_QLMDe}} Class
#' 
#' @description 
#' Quantile least Mahalanobis distance estimates (\code{\link{QLMDe}}) of finite mixture distribution, with the one-dimensional observations included.
#' The \code{\linkS4class{fmx_QLMDe}} object contains (i.e., inherits from) the \code{\linkS4class{fmx}} object. 
#' 
#' @slot data \code{\link[base]{numeric}} vector of the one-dimensional observations
#' 
#' @slot data.name \code{\link[base]{character}} scalar for the name of observations
#' 
#' @slot epdf empirical probability density \code{\link[base]{function}} fitted by \code{\link[stats]{approxfun}}
#' 
#' @slot quantile_vv variance-covariance \code{\link[base]{matrix}} of selected quantiles (based on the selected probabilities stored in slot \code{@@p})
#' 
#' @slot init \code{\linkS4class{fmx}} object indicating the initial values to be sent to \code{\link[stats]{optim}}
#' 
#' @slot p \code{\link[base]{numeric}} vectors of probabilities, where the distance between the empirical and true quantiles will be measured
#' 
#' @slot optim a \code{\link[base]{list}} returned from \code{\link[stats]{optim}}
#' 
#' @examples 
#' ?`fmx_QLMDe-class`
#' 
#' @export
setClass(Class = 'fmx_QLMDe', contains = 'fmx', slots = c(
  data = 'numeric', 
  data.name = 'character',
  epdf = 'function', 
  quantile_vv = 'matrix',
  init = 'fmx',
  p = 'numeric',
  optim = 'list'
), prototype = prototype(
  data.name = 'observations'
), validity = function(object) {
  if (!length(object@data)) stop('*Never* remove @data slot from \'fmx_QLMDe\' object')
  if (anyNA(object@data)) stop('Observations in \'fmx_QLMDe\' must be free of NA_real_')
  if (!length(object@p)) stop('`p` must be recorded')
})



