
#' @title Specification of \linkS4class{fmx} Class
#' 
#' @description 
#' Parameter specification for a one-dimensional finite mixture distribution.
#'
#' @slot distname \link[base]{character} scalar, 
#' name of parametric distribution of the mixture components.
#' Currently, normal (\code{'norm'}) and Tukey's \eqn{g}-&-\eqn{h} (\code{'GH'}) distributions are supported.
#' 
#' @slot parM \link[base]{double} \link[base]{matrix}, 
#' all distribution parameters in the mixture. 
#' Each row corresponds to one component. Each column includes the same parameters of all components.
#' The order of rows corresponds to the (non-strictly) increasing order of the component location parameters.
#' The column names match the formal arguments of the corresponding distribution, 
#' e.g., \code{mean} and \code{sd} for norm distribution (see \code{\link[stats]{dnorm}}), 
#' or \code{A}, \code{B}, \code{g} and \code{h} for Tukey's \eqn{g}-&-\eqn{h} distribution (see \code{\link{dGH}}).
#' 
#' @slot w \link[base]{numeric} vector of mixing proportions that must sum to 1
#' 
# @details ..
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
  if (is.unsorted(parM[,1L], strictly = FALSE)) {
    print(parM[,1L])
    stop('location parameter must be sorted') 
  }
  K <- dim(parM)[1L]
  w <- object@w
  if (K == 1L) {
    if (!isTRUE(all.equal.numeric(w, 1))) stop('slot `w` should be `1` for 1-component') # may not be ?base::identical
  } else {
    if (length(w) != K) stop('slot `w` should be length-K')
    if (!isTRUE(all.equal.numeric(sum(w), 1))) stop('slot `w` must sum up to be `1`') # may not be ?base::identical
  }
  # ?base::is.unsorted.  Note here I use `strictly = FALSE`, 
  # since log(d) may have large negative value, such that exp(log(d)) may be numerically 0.
}) 






