


#' @title A Simpler and Faster Mahalanobis Distance
#' 
#' @description 
#' 
#' Mahalanobis distance, simpler and faster than \code{\link[stats]{mahalanobis}}.
#' 
#' @param x \code{\link[base]{numeric}} vector
#' 
#' @param center \code{\link[base]{numeric}} vector
#' 
#' @param invcov the \emph{inverted} covariance \code{\link[base]{matrix}}
#' 
#' @return 
#' 
#' A \code{\link[base]{numeric}} scalar.
#' 
#' @seealso \code{\link[stats]{mahalanobis}}
#' 
#' @export
mahalanobis_int <- function(x, center, invcov) {
  # if (!is.vector(x, mode = 'double')) stop('x must be double vector') # speed
  # if (!is.vector(center, mode = 'double')) stop('center must be double vector') # speed
  x0 <- x - center
  c(crossprod(x0, invcov) %*% x0)
}


#' @title Inverse of \eqn{X'X} by QR Decomposition
#' 
#' @description Compute \eqn{(X'X)^{-1}} from the R part of the QR decomposition of \eqn{X}.
#' 
#' @param X \eqn{m*n} \code{\link[base]{matrix}}
#' 
#' @seealso \code{\link[base]{chol2inv}}, \code{\link[base]{chol.default}}
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
#' The inverse \code{\link[base]{matrix}} of the cross product \eqn{(X'X)}, where \eqn{X} is the input \code{\link[base]{matrix}}.
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


#' @title Determine Duplicate Elements by \code{\link[base]{all.equal.numeric}}
#' 
#' @description
#' 
#' Determine duplicate elements and extract unique elements in \code{\link[base]{double}} vector, 
#' similar to \code{\link[base]{duplicated}} and \code{\link[base]{unique}},
#' using the algorithm given in \code{\link[base]{all.equal.numeric}}.
#' 
#' @param target \code{\link[base]{double}} vector, missing value \code{NA_real_} is not allowed
#' 
#' @param current \code{\link[base]{double}} vector, missing value \code{NA_real_} is not allowed
#' 
#' @param tolerance see \code{\link[base]{all.equal.numeric}}.
#' 
#' @param x \code{\link[base]{double}} vector
#' 
#' @return 
#' 
#' \code{\link{allequal_dbl}} returns a \code{\link[base]{logical}} \code{\link[base]{matrix}}
#' indicating whether the vector \code{current} is element-wise near-equal to the vector \code{target} 
#' within the prespecified \code{tolerance}.  
#' This function is different from \code{\link[base]{all.equal.numeric}}, such that 
#' (1). only comparisons between real \code{\link[base]{double}} values are performed;
#' (2). element-wise comparison is performed, with the rows of returned \code{\link[base]{matrix}} correspond to \code{current}
#' and columns correspond to \code{target};
#' (3). a \code{\link[base]{logical}} scalar is returned for each element-wise comparison.
#' 
#' \code{\link{duplicated_allequal}} returns a \code{\link[base]{logical}} vector of the same length as the input \code{x},
#' indicating whether each element is a duplicate give the previous elements.  
#' This function behaves otherwise similar to \code{\link[base]{duplicated.default}}, 
#' with the only difference that the comparison is performed by \code{\link{allequal_dbl}}.
#' 
#' \code{\link{unique_allequal}} returns the unique (determined by \code{\link{allequal_dbl}}) elements in the input \code{x},
#' behaving otherwise similar to \code{\link[base]{unique.default}}.
#' 
#' @examples 
#' x = c(.3, 1-.7, 0, .Machine$double.eps)
#' allequal_dbl(current = x, target = c(.3, 0))
#' unique.default(x) # not desired
#' unique_allequal(x) # desired
#' 
#' @name allequal
#' @export
allequal_dbl <- function(target, current, tolerance = sqrt(.Machine$double.eps)) {
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
#' @rdname allequal
#' @export
unique_allequal <- function(x, tolerance = sqrt(.Machine$double.eps)) x[!duplicated_allequal(x, tolerance = tolerance)]
#' @rdname allequal
#' @export
duplicated_allequal <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  if (!is.vector(x, mode = 'double')) stop('input must be double vector')
  nx <- length(x)
  id <- allequal_dbl(target = x, current = x, tolerance = tolerance)
  id[upper.tri(id, diag = TRUE)] <- FALSE # super smart!
  return(.rowSums(id, m = nx, n = nx) > 0) # takes care of 1st element too 
}



