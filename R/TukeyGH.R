

#' @title Tukey's \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description 
#' 
#' Density, distribution function, quantile function and random generation 
#' for the Tukey's \eqn{g}-&-\eqn{h} distribution with 
#' location parameter \eqn{A},
#' scale parameter \eqn{B},
#' skewness \eqn{g} and 
#' kurtosis \eqn{h}.
#' 
#' @param x,q \link[base]{double} \link[base]{vector}, quantiles
#' 
#' @param p \link[base]{double} \link[base]{vector}, probabilities
#' 
#' @param n \link[base]{integer} scalar, number of observations
#' 
#' @param z \link[base]{double} \link[base]{vector}, standard normal quantiles.
#' 
#' @param log,log.p \link[base]{logical} scalar, if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' 
#' @param lower.tail \link[base]{logical} scalar, if \code{TRUE} (default), probabilities are \eqn{Pr(X\le x)} otherwise, \eqn{Pr(X>x)}.
#' 
#' @param A \link[base]{double} scalar, location parameter \eqn{A}, default \eqn{A=0} (as parameter \code{mean} of \link[stats]{dnorm} function)
#' 
#' @param B \link[base]{double} scalar, scale parameter \eqn{B>0}, default \eqn{B=1} (as parameter \code{sd} of \link[stats]{dnorm} function)
#' 
#' @param g \link[base]{double} scalar, skewness parameter \eqn{g}, default \eqn{g=0} indicating no skewness
#' 
#' @param h \link[base]{double} scalar, kurtosis parameter \eqn{h\geq 0}, default \eqn{h=0} indicating no kurtosis
#' 
#' @param q0 \link[base]{double} \link[base]{vector} of \eqn{(q-A)/B}, for internal use to increase compute speed
#' 
# @param interval interval of standard normal quantiles, when solving from Tukey \eqn{g}-&-\eqn{h} quantiles using the vuniroot algorithm 
#' 
#' @param ... other parameters of \link{vuniroot2}
#' 
#' @details
#' 
#' Argument \code{A}, \code{B}, \code{g} and \code{h} will be recycled to the maximum length of the four.
#' 
#' @return 
#' 
#' \link{dGH} gives the density and accommodates \link[base]{vector} arguments \code{A}, \code{B}, \code{g} and \code{h}.
#' The quantiles \code{x} can be either \link[base]{vector} or matrix.
#' This function takes about 1/5 time of \link[gk]{dgh}.
#' 
#' \link{pGH} gives the distribution function, only taking scalar arguments and \link[base]{vector} quantiles \code{q}.
#' This function takes about 1/10 time of \link[gk]{pgh} and \link[OpVaR]{pgh} functions.
#' 
#' \link{qGH} gives the quantile function, only taking scalar arguments and \link[base]{vector} probabilities \code{p}.
#' This function takes about 1/2 time of \link[gk]{qgh} and 1/10 time of \link[OpVaR]{qgh} functions.
#' 
#' \link{rGH} generates random deviates, only taking scalar arguments.
#' 
#' \link{z2qGH} is the Tukey's \eqn{g}-&-\eqn{h} transformation.
#' Note that \code{gk:::z2gh} is only an \strong{approximation} to Tukey's \eqn{g}-&-\eqn{h} transformation.
#' 
#' Unfortunately, \link{qGH2z}, the inverse of Tukey's \eqn{g}-&-\eqn{h} transformation, 
#' does not have a closed form and needs to be solved numerically.
#' 
#' @seealso \link[OpVaR]{dgh} \link[gk]{dgh}
#' 
#' 
#' @examples
#' 
#' (x = c(NA_real_, rGH(n = 5L, g = .3, h = .1)))
#' dGH(x, g = c(0,.1,.2), h = c(.1,.1,.1))
#' 
#' p0 = seq.int(0, 1, by = .2)
#' (q0 = qGH(p0, g = .2, h = .1))
#' range(pGH(q0, g = .2, h = .1) - p0)
#'
#' q = (-2):3; q[2L] = NA_real_; q
#' (p1 = pGH(q, g = .3, h = .1))
#' range(qGH(p1, g = .3, h = .1) - q, na.rm = TRUE)
#' (p2 = pGH(q, g = .2, h = 0))
#' range(qGH(p2, g = .2, h = 0) - q, na.rm = TRUE)
#' 
#' curve(dGH(x, g = .3, h = .1), from = -2.5, to = 3.5)
#' 
#' @name TukeyGH
#' @export
dGH <- function(x, A = 0, B = 1, g = 0, h = 0, log = FALSE, ...) {
  pars <- cbind(A, B, g, h) # recycle
  return(.dGH(x = x, A = pars[,1L], B = pars[,2L], g = pars[,3L], h = pars[,4L], log = log, ...))
}

# not compute intensive
.dGH <- function(x, A, B, g, h, log, interval = c(-50, 50), tol = .Machine$double.eps^.25, maxiter = 1000) {
  # use wider `interval` since not compute intensive
  if (!(nx <- length(x))) return(double(length = 0L)) # ?fitdistrplus::fitdist will test len-0 `x`
  nA <- length(A)
  nB <- length(B)
  ng <- length(g)
  nh <- length(h)

  xok <- is.finite(x) # ?fitdistrplus::fitdist will test exceptions of x = c(0, 1, Inf, NaN, -1)
  
  if ((nA == 1L) && (nB == 1L) && (ng == 1L) && (nh == 1L)) {
    z <- x
    if ((h < 0) || (B < 0)) { # exception handling for ?fitdistrplus::fitdist
      z[] <- NaN
      return(z)
    }
    z[xok] <- .qGH2z(q = c(x[xok]), A = A, B = B, g = g, h = h, interval = interval, tol = tol, maxiter = maxiter)
    
  } else if ((nA == nB) && (nA == ng) && (nA == nh)) {
    #if (!all(xok)) stop('my fmx algorithm do not allow NA or Inf quantile')
    if (is.matrix(x)) {
      if (dim(x)[1L] != nA) stop('nrow of `x` do not match length of `A`')
      z <- q0 <- (x - A)/B
    } else if (is.numeric(x)) {
      z <- q0 <- tcrossprod(1/B, x) - A/B
    } else stop('illegal x: ', sQuote(class(x)[1L]))
    qok <- is.finite(q0) # not `xok` when `x` ?base::is.vector
    for (i in seq_len(nA)) {
      iok <- qok[i,]
      z[i,iok] <- .qGH2z(q0 = q0[i,iok], g = g[i], h = h[i], interval = interval, tol = tol, maxiter = maxiter)
    }
    
  } else stop('length of parameters must match')
  
  if (any(id <- is.infinite(z))) { # `z` is either vector or 'matrix'
    z[id & (z < 0)] <- interval[1L]
    z[id & (z > 0)] <- interval[2L]
  }
  
  ret_log <- -z^2/2 - log(2*pi)/2 - Deriv_z2qGH(z, B = B, g = g, h = h)
  if (log) return(ret_log)
  return(exp(ret_log))
  
}


# Derivative of \link{z2qGH} against `z`, on the log-scale
# inspired by ?OpVaR:::deriv_gh
# Inf in `z` \strong{will} cause trouble
# not sure of the usage of ?base::tanh and ?base::cosh in ?gk:::Qgh_deriv
Deriv_z2qGH <- function(z, B, g, h) {
  hz2 <- h * z^2
  if (length(g) == 1L) { # length(B) == length(h) == 1L; is.vector(z, mode = 'numeric')
    if (g == 0) {
      trm2 <- 1 + hz2
    } else {
      e_gz <- exp(g*z)
      trm2 <- e_gz + h * z * (e_gz - 1)/g
    }
  } else { # length(B) == length(g) == length(h); is.matrix(z); nrow(z) = length(B)
    g1 <- (g != 0)
    z_g1 <- z[g1, , drop = FALSE]
    e_gz1 <- exp(g[g1] * z_g1)
    trm2 <- 1 + hz2 # for `g == 0`, also create 'array'
    trm2[g1,] <- e_gz1 + h[g1] * z_g1 * (e_gz1 - 1)/g[g1]
  }
  
  return(log(B) + hz2/2 + log(trm2))
}







# not compute-intensive
#' @rdname TukeyGH
#' @export
rGH <- function(n, A = 0, B = 1, g = 0, h = 0) z2qGH(rnorm(n), A = A, B = B, g = g, h = h)



#' @rdname TukeyGH
#' @export
qGH <- function(p, A = 0, B = 1, g = 0, h = 0, lower.tail = TRUE, log.p = FALSE) {
  # only works with vector `p` and len-1 `A`,`B`,`g`,`h`, for now
  z <- qnorm(p = p, lower.tail = lower.tail, log.p = log.p)
  .z2qGH(z = z, A = A, B = B, g = g, h = h)
}


# not compute-intensive :)
#' @rdname TukeyGH
#' @export
pGH <- function(q, A = 0, B = 1, g = 0, h = 0, lower.tail = TRUE, log.p = FALSE, ...) {
  # only works with vector `q` and len-1 `A`,`B`,`g`,`h`
  z <- qGH2z(q = q, A = A, B = B, g = g, h = h, ...)
  pnorm(q = z, mean = 0, sd = 1, lower.tail = lower.tail, log.p = log.p)
}





#' @rdname TukeyGH
#' @export
z2qGH <- function(z, A = 0, B = 1, g = 0, h = 0) {
  # not compute intensive
  # z must be numeric vector (i.e. not 'matrix')
  nA <- length(A)
  nB <- length(B)
  ng <- length(g)
  nh <- length(h)
  if (nA == 1L && nB == 1L && ng == 1L && nh == 1L) return(.z2qGH(z, A = A, B = B, g = g, h = h))
  if ((nA != nB) || (nA != ng) || (nA != nh)) stop('distribution parameters must be of same length')
  g0 <- (g == 0) # vector
  h0 <- (h == 0)
  q <- exp(tcrossprod(h, z^2/2))
  q[g0,] <- tcrossprod(g0[g0], z) * q[g0, , drop = FALSE] # \strong{not} \code{g[g0]}
  q[!g0,] <- expm1(tcrossprod(g[!g0], z)) / g[!g0] * q[!g0, , drop = FALSE]
  return(A + q * B)
}


# Tukey GH definition/transformation; not compute intensive
.z2qGH <- function(z, A = 0, B = 1, g = 0, h = 0) {
  # `z` must be numeric vector (i.e. not 'matrix'), all arguments `A`, `B`, `g` and `h` are len-1
  g0 <- (g == 0) # scalar
  h0 <- (h == 0)
  q <- if (g0 && h0) {
    z
  } else if (g0 && !h0) {
    z * exp(h * z^2/2) 
  } else if (!g0 && h0) {
    expm1(g*z) / g
    # if ((h == 0) && (g > 0)) then q > -1/g (when z \arrow -Inf)
    # if ((h == 0) && (g < 0)) then q < -1/g (when z \arrow Inf)
    # numerically, the threshold of epsilon (for `h`) will depend on `g`
  } else {
    expm1(g*z) / g * exp(h * z^2/2)
  }
  return(A + q * B)
}




# not compute intensive (for compute intensive jobs, use \link{.qGH2z})
# inspired by ?OpVaR:::gh_inv
#' @rdname TukeyGH
#' @export
qGH2z <- function(q, q0 = (q - A)/B, A = 0, B = 1, ...) {
  # ?base::is.finite finds finite AND non-missing; as fast as \code{rep(TRUE, times = nq)} (where nq = length(q))
  if (!length(q0)) return(numeric()) # required by ?fitdistrplus::fitdist
  out <- q0
  qok <- is.finite(q0)
  out[qok] <- .qGH2z(q0 = q0[qok], ...)
  return(out)
}




if (FALSE) {
  z = rnorm(1e3L)
  all.equal.numeric(.qGH2z(z2qGH(z, g = .3, h = .1), g = .3, h = .1), z)
  all.equal.numeric(.qGH2z(z2qGH(z, g = 0, h = .1), g = 0, h = .1), z)
  all.equal.numeric(.qGH2z(z2qGH(z, g = .2, h = 0), g = .2, h = 0), z)
}


# internal workhorse of \link{qGH2z}
# inverse of Tukey GH transformation; compute intensive!!!
.qGH2z <- function(
  q, q0 = (q - A)/B, # `q` and `q0` both finite AND non-missing
  A = 0, B = 1, g = 0, h = 0, # all len-1 (`A` and `B` not needed if `q0` is provided)
  interval = c(-15, 15), # smaller `interval` for \code{QLMDe} algorithm (support of standard normal distribution)
  tol = .Machine$double.eps^.25, maxiter = 1000
) {
  
  #if (!length(q0)) return(numeric()) # required by \link[fitdistrplus]{fitdist}
  g0 <- (g == 0)
  h0 <- (h == 0)
  out <- q0
  #interval <- t.default(array(interval, dim = c(2L, length(q0))))
  
  # bound issue only in \code{dGH}, not \code{qfmx}
  
  if (!g0 && !h0) { # most likely to happen in ?stats::optim; put in as the first option to save time
    out[] <- vuniroot2(y = q0 * g, f = function(z) expm1(g*z) * exp(h * z^2/2), interval = interval, tol = tol, maxiter = maxiter)
    # very small `h` would cause bound-issue
    return(out)
  }
  
  if (g0 && h0) return(q0)
  
  if (!g0 && h0) { # has bound but also has explicit form!
    egz <- q0 * g + 1
    if (any(id <- (egz <= 0))) {
      out[id] <- if (g < 0) Inf else -Inf
    }
    out[!id] <- log(egz[!id]) / g
    return(out)
  }
  
  if (g0 && !h0) { # wont have the bound issue if g0
    out[] <- vuniroot2(y = q0, f = function(z) z * exp(h * z^2/2), interval = interval, tol = tol, maxiter = maxiter)
    return(out)
  }
  
}




# need Tingting's \link{ggcurve} function
# ggcurve(dGH, g = .3, h = 0, xlim = c(-3, 5)) # need Tingting's \link{dGH} function
# # learn Hoaglin (1985)
# gv = (1:5)/5
# hv = (1:5)/10
# ggcurve(\(z,g) expm1(g*z)/g, g = gv, xlim = c(-2,2), title = 'Fig 11-2')
# ggcurve(\(z,h) dGH(z,g=0,h=h)-dnorm(z), h = hv, xlim = c(0,8), title = 'Fig 11-8')
# ggcurve(\(z,h) z*exp(h*z^2/2), h = hv, xlim = c(-2,2), title = 'Fig 11-9')
# ggcurve(\(z,h) z*exp(h*z^2/2), h = -hv, xlim = c(-6,6), title = 'Not monotone')






