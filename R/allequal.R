

#' @title Determine Nearly-Equal Elements
#' 
#' @description 
#' Determine nearly-equal elements and extract non-nearly-equal elements in a \link[base]{double} \link[base]{vector}.
#' 
#' @param x \link[base]{double} vector
#' 
#' @param ... additional parameters of function [outer_allequal()]
#' 
#' @returns 
#' 
#' [duplicated_allequal] returns a \link[base]{logical} \link[base]{vector} of the same length as the input vector,
#' indicating whether each element is nearly-equal to any of the previous elements.  
#' 
#' [unique_allequal] returns the non-nearly-equal elements in the input vector.
#' 
#' @seealso 
#' [outer_allequal()] 
#' \link[base]{duplicated.default} \link[base]{unique.default}
#' 
#' @examples 
#' x = c(.3, 1-.7, 0, .Machine$double.eps)
#' unique.default(x) # not desired
#' unique_allequal(x) # desired
#' 
#' @name ud_allequal
#' @export
unique_allequal <- function(x, ...) x[!duplicated_allequal(x, ...)]


#' @rdname ud_allequal
#' @export
duplicated_allequal <- function(x, ...) {
  x <- unclass(x)
  if (!is.vector(x, mode = 'double')) stop('input must be double vector')
  nx <- length(x)
  id <- outer_allequal(target = x, current = x, ...)
  id[upper.tri(id, diag = TRUE)] <- FALSE # super smart!
  return(.rowSums(id, m = nx, n = nx, na.rm = FALSE) > 0) # takes care of 1st element too 
}



