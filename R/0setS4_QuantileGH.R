
#' @title Specification of S4 Object \code{'fmx'}
#' 
#' @description 
#' An S4 class representing a one-dimensional finite mixture distribution.
#'
#' @slot distname 'character' value for the name of parametric distribution of mixture components.
#' Currently normal (\code{'norm'}) and Tukey's \eqn{g}-&-\eqn{h} (\code{'GH'}) distributions are supported.
#' 
#' @slot parM 'matrix' of all distribution parameters in the mixture. 
#' Each row corresponds to one component. Each column includes the same parameters of all components.
#' The order of rows corresponds to the (non-strictly) increasing order of the component location parameters.
#' The column names, e.g., \code{'mean'} and \code{'sd'} for norm distribution, or \code{'A'}, \code{'B'}, \code{'g'} and \code{'h'} for Tukey's \eqn{g}-&-\eqn{h} distribution,
#' are the same as argument names used in \code{\link[stats]{dnorm}} or our \code{\link{dGH}} functions, respectively.
#' 
#' @slot w 'numeric' vector of mixing proportions that must sum to 1
#' 
#' @details ..
#' 
#' @examples 
#' ?`fmx-class`
#' 
#' @export
setClass(Class = 'fmx', slots = c(
  distname = 'character',
  parM = 'matrix',
  w = 'numeric'
), validity = \(object) {
  if (anyNA(object@parM)) stop('do not allow NA in `fmx` distribution parameter')
  if (.Internal(is.unsorted(object@parM[,1L], FALSE))) stop('location parameter must be sorted') 
  # ?base::is.unsorted.  Note here I use `strictly = FALSE`, 
  # since log(d) may have large negative value, such that exp(log(d)) may be numerically 0.
}) 



#' @title Defining S4 Object \code{'fmx_QLMDe'}
#' 
#' @description 
#' Quantile least Mahalanobis distance estimates (\code{\link{QLMDe}}) of finite mixture distribution, with the one-dimensional observations included.
#' The \code{'fmx_QLMDe'} object contains (i.e., inherits from) the \code{'fmx'} object. 
#' 
#' @slot data 'numeric' vector of the one-dimensional observations
#' 
#' @slot data.name 'character' scalar for the name of observations
#' 
#' @slot epdf empirical probability density functions fitted by \code{\link[stats]{approxfun}}
#' 
#' @slot quantile_vv variance-covariance matrix of selected quantiles (based on the selected probabilities \code{@@p})
#' 
#' @slot init \code{'fmx'} object indicating the initial values to be sent to \code{\link[stats]{optim}}
#' 
#' @slot p 'numeric' vectors of probabilities, where the distance between the empirical and true quantiles will be measured
#' 
#' @slot optim 'list' object returned from \code{\link[stats]{optim}}
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
  
  # deprecated slot (no need to preserved for simulation results)
  #constraint = 'character', 'character' vector of \eqn{g} and \eqn{h} parameters to be fixed at zero, e.g., \code{c('g1', 'h2')}.
  #id_constraint = 'logical', 'logical' vector carrying the same information as slot \code{@@constraint}, primarily used by developer
  #type = 'character', # 'mahalanobis' (in use) or 'euclidean' (deprecated)
  #logLik = 'logLik', # log-likelihood, see ?stats::logLik
  #obs_p = 'numeric', # cumulative density of the observations, based on this model fit
  #parTrace = 'matrix', # no longer saved in object@optim$parTrace (but this is S3 and easy to revert)
  #vcov = 'matrix' # deprecated slot
  # end of deprecated slot
  
), prototype = prototype(
  data.name = 'observations'
), validity = \(object) {
  if (!length(object@data)) stop('*Never* remove @data slot from \'fmx_QLMDe\' object')
  if (anyNA(object@data)) stop('Observations in \'fmx_QLMDe\' must be free of NA_real_') # len-0 compatible
  if (!length(object@p)) stop('`p` must be recorded')
})



