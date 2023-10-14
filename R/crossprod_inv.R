


#' @title Inverse of \eqn{X'X} by \eqn{QR} Decomposition
#' 
#' @description 
#' Compute \eqn{(X'X)^{-1}} from the \eqn{R} part of the \eqn{QR} decomposition of \eqn{X}.
#' 
#' @param X \eqn{m*n} \link[base]{matrix} \eqn{X}
#' 
#' @returns 
#' 
#' Function [crossprod_inv()] returns the inverse \link[base]{matrix} of cross product \eqn{X'X}.
#' 
#' @references 
#' \url{https://en.wikipedia.org/wiki/QR_decomposition}, section **Rectangular matrix**
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
#' @seealso 
#' \link[base]{chol2inv} \link[base]{chol.default}
#' \link[base]{qr.default} \link[base]{qr.R} \link[base]{chol2inv}
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



