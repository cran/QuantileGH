

# Note that we do not have parameter `...` in [vuniroot2], otherwise it will cause warning message 
# `The [dGH]/[pGH] function should raise an error when names are incorrectly named' in `fitdistrplus:::test1fun`.


#' @title Vectorised One Dimensional Root (Zero) Finding
#' 
#' @description 
#' 
#' To solve a monotone function \eqn{y = f(x)} for a given \link[base]{vector} of \eqn{y} values.
#' 
#' @param y \link[base]{numeric} \link[base]{vector} of \eqn{y} values
#' 
#' @param f monotone \link[base]{function} \eqn{f(x)} whose roots are to be solved
#' 
#' @param interval length two \link[base]{numeric} \link[base]{vector}
#' 
#' @param tol \link[base]{double} scalar, desired accuracy (convergence tolerance), 

#' @param maxiter \link[base]{integer} scalar, maximum number of iterations
#' 
#' @details
#' 
#' Function [vuniroot2()], different from \link[rstpm2]{vuniroot}, does
#' \itemize{
#' \item{accept `NA_real_` as element(s) of \eqn{y}}
#' \item{handle the case when the analytic root is at `lower` and/or `upper`}
#' \item{return a root of 
#' `Inf` (if `abs(f(lower)) >= abs(f(upper))`) or 
#' `-Inf` (if `abs(f(lower)) < abs(f(upper))`), 
#' when the function value `f(lower)` and `f(upper)` are not of opposite sign.}
#' }
#' 
#' @returns 
#' 
#' Function [vuniroot2()] returns a \link[base]{numeric} \link[base]{vector} \eqn{x} as the solution of \eqn{y = f(x)} with given vector \eqn{y}.
#' 
#' @examples 
#' library(rstpm2)
#' lwr = rep(1, times = 9L); upr = rep(3, times = 9L)
#' 
#' # ?rstpm2::vuniroot does not accept NA \eqn{y}
#' tryCatch(vuniroot(function(x) x^2 - c(NA, 1:8), lower = lwr, upper = upr), error = identity)
#' 
#' # ?rstpm2::vuniroot not good when the analytic root is at `lower` or `upper`
#' f <- function(x) x^2 - 1:9
#' tryCatch(vuniroot(f, lower = lwr, upper = upr, extendInt = 'no'), warning = identity)
#' tryCatch(vuniroot(f, lower = lwr, upper = upr, extendInt = 'yes'), warning = identity)
#' tryCatch(vuniroot(f, lower = lwr, upper = upr, extendInt = 'downX'), error = identity)
#' tryCatch(vuniroot(f, lower = lwr, upper = upr, extendInt = 'upX'), warning = identity)
#' 
#' vuniroot2(c(NA, 1:9), f = function(x) x^2, interval = c(1, 3)) # all good
#' 
#' @importFrom rstpm2 vuniroot
#' @export
vuniroot2 <- function(
  y, f, interval = stop('must provide a length-2 `interval`'), 
  tol = .Machine$double.eps^.25, maxiter = 1000L
) {
  
  if (any(is.infinite(y))) stop('infinite return from function `f` cannot be handled')
  out <- rep(NA_real_, times = length(y)) # y * NA_real_; # since NA * Inf -> NaN
  if (!any(yok <- !is.na(y))) return(out)
  out[yok] <- Inf
  yok_ <- y[yok]
  
  f.intv <- f(interval) # len-2
  if (anyNA(f.intv)) {
    #f <<- f; interval <<- interval
    stop('function evaluated at either end of the `interval` must not be NA')
  }
  #if (any(is.infinite(f.intv))) { ### Inf ends are fine!!
  #  #f <<- f; interval <<- interval
  #  #stop('function evaluated at either end of the `interval` must not be Inf or -Inf')
  #  return(out)
  #}
  
  # Remove this.  Just let algorithm declare the solution being Inf
  #if (f.intv[1L] == f.intv[2L]) stop('function in `vuniroot2` not monotone')
  
  f.lower <- f.intv[1L] - yok_
  f.upper <- f.intv[2L] - yok_
  
  # check values at either end of `interval`
  if (any(fl0 <- (abs(f.lower) < tol))) {
    out[yok][fl0] <- interval[1L]
    f.lower <- f.lower[!fl0]
    f.upper <- f.upper[!fl0]
    yok[yok][fl0] <- FALSE # update `yok` the last
    if (!any(yok)) return(out)
  }
  if (any(fu0 <- (abs(f.upper) < tol))) {
    out[yok][fu0] <- interval[2L]
    f.lower <- f.lower[!fu0]
    f.upper <- f.upper[!fu0]
    yok[yok][fu0] <- FALSE # update `yok` the last
    if (!any(yok)) return(out)
  }
  yok_ <- y[yok]
  # end of checking at either end of `interval`

  if (any(sign_same <- (f.lower * f.upper > 0))) {
    # smart!  used in ?base::uniroot and ?rstpm2::vuniroot
    # `sign_same` are the indices of `y` where {f(interval[1L]) - y} and {f(interval[2L]) - y} are of the same signs
    # they will cause error in ?rstpm2::vuniroot and/or ?base::uniroot
    # therefore I will let the solution at these indices to be either -Inf or Inf
    out[yok][sign_same & (abs(f.lower) < abs(f.upper))] <- -Inf # otherwise Inf (as defined as default)
    sign_change <- which(!sign_same)
    # `sign_change` are the indices of `y` where {f(interval[1L]) - y} and {f(interval[2L]) - y} are of opposite signs
  } else sign_change <- seq_along(yok_)
  
  if (nn <- length(sign_change)) { 
    # very important!! 
    # `sign_change` could be len-0
    # in such case R will crash (instead of throw an ?rstpm2::vuniroot error)
    # ?base::tryCatch will not be able to stop R from crashing!
    out[yok][sign_change] <- vuniroot(
      f = function(x) f(x) - yok_[sign_change],
      lower = rep(interval[1L], times = nn), upper = rep(interval[2L], times = nn), 
      f.lower = f.lower[sign_change], f.upper = f.upper[sign_change], 
      extendInt = 'no', check.conv = TRUE,
      tol = tol, maxiter = maxiter, trace = 0L)[[1L]] # do *not* ?base::suppressWarnings
  }
  
  return(out)
  
}

