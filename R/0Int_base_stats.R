


# @title Console output in color
# @description Inspired by \code{\link[insight]{color_text}}
# @param x 'character' vector, 0-char accepted, while \code{NA_character_} should be removed (won't give error though)
# @seealso \code{insight:::.colour}
# @examples
# cat(.cyan(c('hello', 'world')), '\n')
# cat(.bold(.red(c('hello', 'world'))), '\n')
# @name insight_color
# @export
.cyan <- function(x) paste0('\033[36m', x, '\033[39m') # ?insight:::.cyan
# @rdname insight_color
# @export
.red <- function(x) paste0('\033[31m', x, '\033[39m') # ?insight:::.red
# @rdname insight_color
# @export
.yellow <- function(x) paste0('\033[33m', x, '\033[39m') # ?insight:::.yellow
# @rdname insight_color
# @export
.green <- function(x) paste0('\033[32m', x, '\033[39m') # ?insight:::.green
# @rdname insight_color
# @export
.blue <- function(x) paste0('\033[34m', x, '\033[39m') # ?insight:::.blue
# @rdname insight_color
# @export
.violet <- function(x) paste0('\033[35m', x, '\033[39m') # ?insight:::.violet
# @rdname insight_color
# @export
.grey <- function(x) paste0('\033[90m', x, '\033[39m') # ?insight:::.grey
# @rdname insight_color
# @export
.bold <- function(x) paste0('\033[1m', x, '\033[22m') # ?insight:::.bold
# @rdname insight_color
# @export
.italic <- function(x) paste0('\033[3m', x, '\033[23m') # ?insight:::.italic








# speeds up ?stats::mahalanobis
mahalanobis_int <- function(x, center, invcov) {
  # if (!is.vector(x, mode = 'double')) stop('x must be double vector') # speed
  # if (!is.vector(center, mode = 'double')) stop('center must be double vector') # speed
  # `invcov` is the \strong{inverted} covariance matrix
  x0 <- x - center
  #c(.Internal(crossprod(x0, invcov)) %*% x0) # .Internal needed for speed
  t.default(x0) %*% invcov %*% x0 # not that slow in practice.. remove .Internal for CRAN
}


# based on z- or t-test; see ?stats::confint.default
confint_int <- function(cf, vv, ses, df = Inf, level = .95, ...) {
  if (missing(ses)) {
    if (isS4(vv)) vv <- as.matrix(vv)
    if (!is.matrix(vv)) stop('`vv` must be convertible to \'matrix\'')
    if ((dmv <- dim(vv))[1L] != dmv[2L]) stop('`vv` needs to be square matrix')
    ses <- sqrt(vv[cbind(vseq <- seq_len(dmv[1L]), vseq)]) # \strong{not} matrix squared-root!
  } # else { dont bother to check length(ses) == length(cf)}
  a <- (1 - level)/2 # stopifnot(length(level) == 1L) # save time
  ndf <- length(df)
  fac <- qt(a, df = df) # is.infinite(df) compatible
  ci <- if (ndf == 1L) {
    cf + tcrossprod(ses, c(fac, -fac)) # remove .Internal for CRAN
  } else if (ndf == length(ses)) {
    cf + ses * cbind(fac, -fac)
  } else stop('length(df) wrong')
  dimnames(ci) <- list(names(cf), sprintf(fmt = '%.1f%%', 1e2*c(a, 1-a)))
  ci
}



#' @title Inverse of \eqn{X'X} by QR Decomposition
#' 
#' @description Compute \eqn{(X'X)^{-1}} from the R part of the QR decomposition of \eqn{X}.
#' 
#' @param X \eqn{m*n} matrix
#' 
#' @seealso \code{\link[base]{chol2inv}}
#' 
#' @references 
#' \url{https://en.wikipedia.org/wiki/QR_decomposition}, section \strong{Rectangular matrix}
#' 
#' @examples 
#' set.seed(123); (X = array(rnorm(40L), dim = c(8L, 5L)))
#' (Xcp = crossprod(X))
#' all.equal.numeric(solve(Xcp), crossprod_inv(X))
#' 
# \donttest{
# library(microbenchmark)
# microbenchmark(solve.default(Xcp), crossprod_inv(X))
# for (i in 1:1e5L) {
#  X = array(rnorm(40L), dim = c(8L, 5L))
#  stopifnot(all.equal.numeric(solve(crossprod(X)), crossprod_inv(X)))
# }
# }
#' 
#' @return 
#' 
#' The inverse \code{'matrix'} of the cross product \eqn{(X'X)}, where \eqn{X} is the input \code{'matrix'}.
#' 
#' @export
crossprod_inv <- function(X) {
  Xqr <- qr.default(X) # ?base::qr.default uses Fortran code
  XR <- qr.R(Xqr) # ?base::qr.R could be slow too
  # all.equal.numeric(current = qr.Q(Xqr) %*% XR, target = X)
  # all.equal.numeric(target = .Internal(crossprod(X, X)), current = .Internal(crossprod(XR, XR)))
  # `chol(.Internal(crossprod(X, X)))` and `XR` not the same (?base::abs same)
  ret <- chol2inv(XR)
  # all.equal.numeric(current = ret, target = chol2inv(chol(.Internal(crossprod(XR, XR)))))
  return(ret)
}


# @title Determine Duplicate Elements by \code{\link[base]{all.equal.numeric}}
# @description ..
# @param target,current,tolerance see \code{\link[base]{all.equal.numeric}}
# @param x 'numeric' vector
# @param ... ..
# @examples 
# x = c(.3, 1-.7, 0, .Machine$double.eps)
# unique.default(x)
# unique_allequal(x)
# \dontrun{
# library(microbenchmark)
# x1 = rep(x, times = 10L)
# microbenchmark(duplicated.default(x), duplicated_allequal(x))
# }
# @return 
# \code{\link{equal_dbl}} returns .. 
# \code{\link{unique_allequal}} returns ..
# \code{\link{duplicated_allequal}} returns ..
# @name equal_dbl
# @export
equal_dbl <- function(target, current, tolerance = sqrt(.Machine$double.eps), ...) {
  # `target` and `current` both 'double'
  # simplified from ?base::all.equal.numeric
  if (length(target) != 1L || !(nc <- length(current))) stop('target must be len-1')
  tok <- !is.na(target)
  cok <- !is.na(current)
  if (!tok) return(!cok) # mimic \code{duplicated.default(c(NA, NA))}
  
  out <- logical(length = nc)
  xy <- abs(target - current[cok]) # 'absolute'
  xn <- abs(target)
  if (is.finite(xn) && xn > tolerance) xy <- xy/xn # 'relative'
  out[cok][xy <= tolerance] <- TRUE
  return(out)
}


# @rdname equal_dbl
# @export
unique_allequal <- function(x, ...) x[!duplicated_allequal(x, ...)]


# @rdname equal_dbl
# @export
duplicated_allequal <- function(x, ...) {
  if (!is.vector(x, mode = 'double')) stop('input must be double vector')
  out <- x
  storage.mode(out) <- 'logical'
  out[1L] <- FALSE
  id <- seq_len(length(x))[-1L]
  out[id] <- vapply(id, FUN = \(i) {
    any(equal_dbl(target = x[i], current = x[seq_len(i-1L)], ...))
  }, FUN.VALUE = NA)
  return(out)
}





