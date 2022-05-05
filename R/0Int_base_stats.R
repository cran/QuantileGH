


#' @title A Simpler and Faster Mahalanobis Distance
#' 
#' @description 
#' 
#' A simpler and faster \link[stats:mahalanobis]{Mahalanobis} distance.
#' 
#' @param x \link[base]{numeric} vector
#' 
#' @param center \link[base]{numeric} vector, mean \eqn{\mathbf{\mu}}
#' 
#' @param invcov \link[base]{numeric} \link[base]{matrix}, \emph{inverted} variance-covariance \eqn{\mathbf{\Sigma}}
#' 
#' @return 
#' 
#' \link{mahalanobis_int} returns a \link[base]{numeric} scalar.
#' 
#' @seealso \link[stats]{mahalanobis}
#' 
#' @export
mahalanobis_int <- function(x, center, invcov) {
  # if (!is.vector(x, mode = 'double')) stop('x must be double vector') # speed
  # if (!is.vector(center, mode = 'double')) stop('center must be double vector') # speed
  x0 <- x - center
  c(crossprod(x0, invcov) %*% x0)
}


#' @title Inverse of \eqn{X'X} by \eqn{QR} Decomposition
#' 
#' @description Compute \eqn{(X'X)^{-1}} from the \eqn{R} part of the \eqn{QR} decomposition of \eqn{X}.
#' 
#' @param X \eqn{m*n} \link[base]{matrix}
#' 
#' @seealso \link[base]{chol2inv} \link[base]{chol.default}
#' 
#' @references 
#' \url{https://en.wikipedia.org/wiki/QR_decomposition}, section \strong{Rectangular matrix}
#' 
#' @examples 
#' set.seed(123); (X = array(rnorm(40L), dim = c(8L, 5L)))
#' stopifnot(all.equal.numeric(solve(crossprod(X)), crossprod_inv(X)))
#' 
# \donttest{
# library(microbenchmark)
# microbenchmark(solve.default(crossprod(X)), crossprod_inv(X))
# for (i in 1:1e5L) {
#  X = array(rnorm(40L), dim = c(8L, 5L))
#  stopifnot(all.equal.numeric(solve(crossprod(X)), crossprod_inv(X)))
# }
# }
#' 
#' @return 
#' 
#' \link{crossprod_inv} returns the inverse \link[base]{matrix} of cross product \eqn{X'X}.
#' 
#' @export
crossprod_inv <- function(X) {
  Xqr <- qr.default(X) # ?base::qr.default uses Fortran code
  XR <- qr.R(Xqr) # ?base::qr.R could be slow too
  # all.equal.numeric(current = qr.Q(Xqr) %*% XR, target = X)
  # all.equal.numeric(target = crossprod(X, X), current = crossprod(XR, XR))
  # `chol(crossprod(X, X))` and `XR` not the same (?base::abs same, why?)
  ret <- chol2inv(XR)
  # all.equal.numeric(current = ret, target = chol2inv(chol.default(crossprod(XR, XR))))
  return(ret)
}


#' @title Test if Two \link[base]{double} Vectors are Element-Wise (Nearly) Equal
#' 
#' @description
#' 
#' Test if two \link[base]{double} vectors are element-wise (nearly) equal.
#' 
#' @param current length-\eqn{n} \link[base]{double} vector, the value to be compared with \code{target}, missing value not allowed
#' 
#' @param target length-\eqn{n} \link[base]{double} vector, the target value, missing value not allowed
#' 
#' @param tolerance \link[base]{double} scalar, see \link[base]{all.equal.numeric}.
#' 
#' @param ... potential parameters, currently not in use
#' 
#' @details 
#' 
#' \link{outer_allequal} is different from \link[base]{all.equal.numeric}, such that 
#' (1). only comparisons between real \link[base]{double} values are performed;
#' (2). element-wise comparison is performed, with the rows of returned \link[base]{matrix} correspond to \code{current}
#' and columns correspond to \code{target};
#' (3). a \link[base]{logical} scalar is returned for each element-wise comparison.
#' 
#' @return 
#' 
#' \link{outer_allequal} returns an \eqn{m*n} \link[base]{logical} \link[base]{matrix}
#' indicating whether the length-\eqn{n} vector \code{current} is element-wise near-equal to the length-\eqn{m} vector \code{target} 
#' within the prespecified \code{tolerance}.  
#' 
#' @seealso \link[base]{all.equal.numeric} \link[base]{outer}
#' 
#' @examples 
#' x = c(.3, 1-.7, 0, .Machine$double.eps)
#' outer_allequal(current = x, target = c(.3, 0))
#' 
#' @export
outer_allequal <- function(target, current, tolerance = sqrt(.Machine$double.eps), ...) {
  if (!(nt <- length(target))) stop('len-0 `target`')
  if (!(nc <- length(current))) stop('len-0 `current`')
  if (anyNA(target) || anyNA(current)) stop('Do not allow missingness in `target` or `current`')
  
  mc <- array(current, dim = c(nc, nt))
  mt <- t.default(array(target, dim = c(nt, nc)))
  xy <- abs(mt - mc) # 'absolute'
  xn <- abs(mt)
  if (all(is.finite(xn)) && all(xn > tolerance)) xy <- xy/xn # 'relative'
  return(xy <= tolerance)
}




#' @title Determine Nearly-Equal Elements
#' 
#' @description 
#' Determine nearly-equal elements and extract non-nearly-equal elements in a \link[base]{double} vector.
#' 
#' @param x \link[base]{double} vector
#' 
#' @param ... additional parameters of \link{outer_allequal}
#' 
#' @return 
#' 
#' \link{duplicated_allequal} returns a \link[base]{logical} vector of the same length as the input vector,
#' indicating whether each element is nearly-equal to any of the previous elements.  
#' 
#' \link{unique_allequal} returns the non-nearly-equal elements in the input vector.
#' 
#' @seealso \link[base]{duplicated.default} \link[base]{unique.default}
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
  if (!is.vector(x, mode = 'double')) stop('input must be double vector')
  nx <- length(x)
  id <- outer_allequal(target = x, current = x, ...)
  id[upper.tri(id, diag = TRUE)] <- FALSE # super smart!
  return(.rowSums(id, m = nx, n = nx, na.rm = FALSE) > 0) # takes care of 1st element too 
}



