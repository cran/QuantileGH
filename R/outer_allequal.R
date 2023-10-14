
#' @title Test if Two \link[base]{double} Vectors are Element-Wise (Nearly) Equal
#' 
#' @description
#' 
#' Test if two \link[base]{double} \link[base]{vector}s are element-wise (nearly) equal.
#' 
#' @param current length-\eqn{n_c} \link[base]{double} \link[base]{vector}, 
#' the value(s) to be compared with `target`, missing value not allowed
#' 
#' @param target length-\eqn{n_t} \link[base]{double} \link[base]{vector}, 
#' the target value(s), missing value not allowed
#' 
#' @param tolerance positive \link[base]{double} scalar, default `sqrt(.Machine$double.eps)` 
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details 
#' 
#' Function [outer_allequal()] is different from \link[base]{all.equal.numeric}, such that 
#' \itemize{
#' \item{only compares between \link[base]{double}, not \link[base]{complex}, values}
#' \item{element-wise comparison is performed}
#' \item{a \link[base]{logical} scalar is always returned for each element-wise comparison.}
#' }
#' 
#' @returns 
#' 
#' Function [outer_allequal()] returns an \eqn{n_c\times n_t} \link[base]{logical} \link[base]{matrix}
#' indicating whether the length-\eqn{n_c} \link[base]{vector} `current` 
#' is element-wise near-equal to 
#' the length-\eqn{n_t} \link[base]{vector} `target` 
#' within the pre-specified `tolerance`.  
#' 
#' @seealso 
#' \link[base]{all.equal.numeric} \link[base]{outer}
#' 
#' @examples 
#' x = c(.3, 1-.7, 0, .Machine$double.eps)
#' outer_allequal(current = x, target = c(.3, 0))
#' 
#' @export
outer_allequal <- function(target, current, tolerance = sqrt(.Machine$double.eps), ...) {
  if (!(nt <- length(target)) || !is.double(target)) stop('len-0 `target`')
  if (!(nc <- length(current)) || !is.double(current)) stop('len-0 `current`')
  if (anyNA(target) || anyNA(current)) stop('Do not allow missingness in `target` or `current`')
  
  mc <- array(current, dim = c(nc, nt))
  mt <- t.default(array(target, dim = c(nt, nc)))
  xy <- abs(mt - mc) # 'absolute'
  xn <- abs(mt)
  if (all(is.finite(xn)) && all(xn > tolerance)) xy <- xy/xn # 'relative'
  return(xy <= tolerance)
}
