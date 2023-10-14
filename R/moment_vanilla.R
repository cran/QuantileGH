

#' @title Raw, Central and Standardized Moments, and other Distribution Characteristics
#' 
#' @description
#' Up to 4th order of raw \eqn{E(Y^n)}, central \eqn{E[(Y-\mu)^n]} and 
#' standardized moments \eqn{E[(Y-\mu)^n]/\sigma^n} of the random variable
#' \eqn{Y = (X - \mathtt{location})/\mathtt{scale}}, 
#' as well as the distribution characteristics 
#' (e.g., mean, standard deviation, skewness and excess kurtosis) of the random variable \eqn{X}.
#' 
#' @slot distname \link[base]{character} scalar, name of distribution,
#' e.g., `'norm'` for normal, `'sn'` for skew-normal, `'st'` for skew-\eqn{t}, 
#' and `GH` for Tukey's \eqn{g}-and-\eqn{h} distribution,
#' following the nomenclature of \link[stats]{dnorm}, \link[sn]{dsn}, \link[sn]{dst} and [dGH()]
#' 
#' @slot location,scale \link[base]{numeric} \link[base]{vector}s or scalars, 
#' location and scale parameters
#' 
#' @slot mu \link[base]{numeric} \link[base]{vector} or scalar, 
#' 1st order *raw* moment \eqn{\mu = E(Y)}. 
#' Note that the 1st order central moment \eqn{E(Y-\mu)} and
#' standardized moment \eqn{E(Y-\mu)/\sigma} are 0.
#' 
#' @slot raw2,raw3,raw4 \link[base]{numeric} \link[base]{vector}s or scalars, 
#' 2nd or higher order *raw* moments \eqn{E(Y^n)}, \eqn{n\geq 2}
#' 
#' @slot central2,central3,central4 \link[base]{numeric} \link[base]{vector}s or scalars, 
#' 2nd or higher order *central* moments, \eqn{\sigma^2 = E[(Y-\mu)^2]} and 
#' \eqn{E[(Y-\mu)^n]}, \eqn{n\geq 3}
#' 
#' @slot standardized3,standardized4 \link[base]{numeric} \link[base]{vector}s or scalars, 
#' 3rd or higher order *standardized* moments, 
#' skewness \eqn{E[(Y-\mu)^3]/\sigma^3} and
#' kurtosis \eqn{E[(Y-\mu)^4]/\sigma^4}.
# \eqn{E[(Y-\mu)^n]/\sigma^n}, \eqn{n=5,\cdots}. 
#' Note that the 2nd standardized moment is 1
#' 
#' @slot mean,sd,skewness,kurtosis \link[base]{numeric} \link[base]{vector}s or scalars, 
#' distribution characteristics of random variable \eqn{X}, such as mean, standard deviation, skewness, and excess kurtosis
#' 
#' @details
#' 
#' For \eqn{Y = (X - \mathtt{location})/\mathtt{scale}}, 
#' let \eqn{\mu = E(Y)}, then the second to fourth central moments of \eqn{Y} are,
#' \deqn{E[(Y-\mu)^2] = E(Y^2) - 2\mu E(Y) + \mu^2 = E(Y^2) - \mu^2}
#' \deqn{E[(Y-\mu)^3] = E(Y^3) - 3\mu E(Y^2) + 3\mu^2 E(Y) - \mu^3 = E(Y^3) - 3\mu E(Y^2) + 2\mu^3}
#' \deqn{E[(Y-\mu)^4] = E(Y^4) - 4\mu E(Y^3) + 6\mu^2 E(Y^2) - 4\mu^3 E(Y) + \mu^4 = E(Y^4) - 4\mu E(Y^3) + 6\mu^2 E(Y^2) - 3\mu^4}
#' 
#' The distribution characteristics of \eqn{Y} are,
#' \deqn{\mu_Y = \mu}
#' \deqn{\sigma_Y = \sqrt{E[(Y-\mu)^2]}}
#' \deqn{\mathtt{skewness}_Y = E[(Y-\mu)^3] / \sigma^3_Y}
#' \deqn{\mathtt{kurtosis}_Y = E[(Y-\mu)^4] / \sigma^4_Y - 3}
#' 
#' The distribution characteristics of \eqn{X} are
#' \eqn{\mu_X = \mathtt{location} + \mathtt{scale}\cdot \mu_Y},
#' \eqn{\sigma_X = \mathtt{scale}\cdot \sigma_Y},
#' \eqn{\mathtt{skewness}_X = \mathtt{skewness}_Y}, and
#' \eqn{\mathtt{kurtosis}_X = \mathtt{kurtosis}_Y}.
#' 
#' 
#' @returns
#' Functions [moment_()], [moment_GH()], [moment_sn()], [moment_st()], [moment_norm()] all return a \linkS4class{moment} object.
#' 
#' @note
#' Potential name clash with `e1071::moment`.
#' 
#' @references
#' \url{https://en.wikipedia.org/wiki/Binomial_theorem}
#' \url{https://en.wikipedia.org/wiki/Central_moment}
#' \url{https://en.wikipedia.org/wiki/Standardized_moment}
#' \url{https://en.wikipedia.org/wiki/Skewness}
#' \url{https://en.wikipedia.org/wiki/Kurtosis}
#' 
#' @examples
#' library(ggplot2)
#' 
#' moment_(dist = 'norm', mean = 1.234, sd = .58)
#' 
#' \dontrun{ # requires Tingting's \pkg{QuantileGH}
#' A = 3; B = 1.5; g = .7; h = .1
#' moment_(dist = 'GH', A = A, B = B, g = 0, h = h)
#' moment_(dist = 'GH', A = A, B = B, g = g, h = 0)
#' moment_(dist = 'GH', A = A, B = B, g = g, h = h)}
#' 
#' xi = 2; omega = 1.3; alpha = 3; nu = 6
#' 
#' ggplot() + geom_function(fun = sn::dsn, args = list(
#'   xi = xi, omega = omega, alpha = alpha
#' ), xlim = c(0, 6))
#' moment_(dist = 'sn', xi, omega, alpha)
#' 
#' ggplot() + geom_function(fun = sn::dst, args = list(
#'   xi = xi, omega = omega, alpha = alpha, nu = nu
#' ), xlim = c(0, 6))
#' moment_(dist = 'st', xi, omega, alpha, nu)
#' 
#' 
#' @name moment_
#' @aliases moment-class
#' @export
setClass(Class = 'moment', slots = c(
  location = 'numeric', scale = 'numeric',
  distname = 'character',
  mu = 'numeric',
  raw2 = 'numeric', raw3 = 'numeric', raw4 = 'numeric',
  central2 = 'numeric', central3 = 'numeric', central4 = 'numeric',
  standardized3 = 'numeric', standardized4 = 'numeric',
  mean = 'numeric', sd = 'numeric', skewness = 'numeric', kurtosis = 'numeric'
))







# converts raw moments to central (and standardized) moments, using binomial theorem
moment_int <- function(
    distname, location, scale,
    mu, raw2, raw3, raw4,
    ...
) {
  central2 <- raw2 - mu^2
  sigma <- sqrt(central2)
  central3 <- raw3 - 3*raw2*mu + 2*mu^3
  central4 <- raw4 - 4*raw3*mu + 6*raw2*mu^2 - 3*mu^4
  
  standardized3 <- central3 / sigma^3
  standardized4 <- central4 / sigma^4
  
  new(Class = 'moment', 
      distname = distname, location = location, scale = scale,
      mu = mu, raw2 = raw2, raw3 = raw3, raw4 = raw4,
      central2 = central2, central3 = central3, central4 = central4,
      standardized3 = standardized3, standardized4 = standardized4,
      mean = location + scale * mu, 
      sd = scale * sigma, 
      skewness = standardized3,
      kurtosis = standardized4 - 3)
  
}





#' @rdname moment_
#' 
#' @param A,B,g,h \link[base]{numeric} \link[base]{vector}s or scalars,
#' parameters of Tukey's \eqn{gh} distribution [dGH()]
#' 
#' @references 
#' Raw moments of Tukey's GH: \doi{10.1002/9781118150702.ch11}
#' 
#' @export
moment_GH <- function(A, B, g, h) {
  tmp <- data.frame(A, B, g, h) # recycling
  A <- tmp[[1L]]
  B <- tmp[[2L]]
  g <- tmp[[3L]]
  h <- tmp[[4L]]
  
  g0 <- (g == 0)
  
  # 1st-4th order raw moment E(Y^n), when `g = 0`
  mu <- r3 <- rep(0, times = length(A))
  r2 <- 1 / (1-2*h) ^ (3/2) # (45a), page 502
  r4 <- 3 / (1-4*h) ^ (5/2) # (45b), page 502
  
  if (any(!g0)) {
    # {r}aw moment E(Y^n), when `g != 0`
    r_g <- function(n) { # equation (47), page 503
      tmp <- lapply(0:n, FUN = function(i) {
        (-1)^i * choose(n,i) * exp((n-i)^2 * g^2 / 2 / (1-n*h))
      })
      suppressWarnings(Reduce(f = `+`, tmp) / g^n / sqrt(1-n*h)) # warnings for `h > 1/n`
    }
    
    mu[!g0] <- r_g(1L)[!g0]
    r2[!g0] <- r_g(2L)[!g0]
    r3[!g0] <- r_g(3L)[!g0]
    r4[!g0] <- r_g(4L)[!g0]
  }
  
  moment_int(distname = 'GH', location = A, scale = B, mu = mu, raw2 = r2, raw3 = r3, raw4 = r4)
  
}


#' @rdname moment_
#' 
#' @param xi,omega,alpha \link[base]{numeric} \link[base]{vector}s or scalars, 
#' location, scale and slant parameters for skew-normal \link[sn]{dsn} and 
#' skew-\eqn{t} \link[sn]{dst} distributions
#' 
#' @references
#' Raw moments of skew-normal: \url{https://en.wikipedia.org/wiki/Skew_normal_distribution}
#' 
#' @export
moment_sn <- function(xi, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(2/pi)
  mu <- b * delta
  moment_int(distname = 'sn', location = xi, scale = omega,
             mu = mu,
             raw2 = 1,
             raw3 = - pi/2 *mu^3 + 3*mu,
             raw4 = 3)
}


#' @rdname moment_
#' 
#' @param nu positive \link[base]{numeric} \link[base]{vector} or scalar, 
#' degrees of freedom(s) of skew-\eqn{t} \link[sn]{dst} distribution
#' 
#' @references
#' Raw moments of skew-\eqn{t}: \url{https://arxiv.org/abs/0911.2342}
#' 
#' @export
moment_st <- function(xi, omega, alpha, nu) {
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(nu/pi) * gamma(nu/2 - 1/2) / gamma(nu/2) # equation (29); https://arxiv.org/pdf/0911.2342.pdf
  mu <- b * delta
  moment_int(distname = 'st', location = xi, scale = omega,
             mu = mu,
             raw2 = nu/(nu-2),
             raw3 = mu * (3-delta^2) * nu/(nu-3),
             raw4 = 3*nu^2/(nu-2)/(nu-4))
}


#' @rdname moment_
#' 
#' @param mean,sd \link[base]{numeric} \link[base]{vector}s or scalars, mean and standard deviation 
#' parameters for normal \link[stats]{dnorm} distribution
#' 
#' @references
#' Raw moments of normal: \url{https://en.wikipedia.org/wiki/Normal_distribution} (replace with \eqn{\mu = 0} and \eqn{\sigma = 1})
#' 
#' @export
moment_norm <- function(mean, sd) {
  moment_int(distname = 'norm', location = mean, scale = sd, mu = 0, raw2 = 1, raw3 = 0, raw4 = 3)
}


# stats::dt
# Raw moments of \eqn{t}-distribution: \url{https://en.wikipedia.org/wiki/Student%27s_t-distribution}
# Raw moments of non-central \eqn{t}-distribution: \url{https://en.wikipedia.org/wiki/Noncentral_t-distribution}











#' @rdname moment_
#' 
#' @param dist see **Usage**
#' 
#' @param ... distribution parameters as described in **Arguments** for [moment_.character()], 
#' or place holder for S3 method dispatch for [moment_.fmx()]
#' 
#' @export
moment_ <- function(dist, ...) UseMethod('moment_')




#' @rdname moment_
#' 
#' @details
#' The S3 method dispatch [moment_.character()] obtains the moments and distribution characteristics from
#' the distribution name `dist` and parameters given in `...`.
#' 
#' @export
moment_.character <- function(dist = c('norm', 'GH', 'sn', 'st'), ...) {
  do.call(what = paste0('moment_', match.arg(dist)), args = list(...))
}

#' @rdname moment_
#' 
#' @details
#' The S3 method dispatch [moment_.fmx()] obtains the moments and distribution characteristics of each mixture component of 
#' an \linkS4class{fmx} object.
#' 
#' @export
moment_.fmx <- function(dist, ...) {
  pars <- dist@pars 
  par_nm <- colnames(pars)
  x <- lapply(seq_len(dim(pars)[2L]), FUN = function(i) pars[,i])
  names(x) <- par_nm
  do.call(what = paste0('moment_', dist@distname), args = x)
}



#' @title Show \linkS4class{moment}
#' 
#' @description ..
#' 
#' @param object \linkS4class{moment}
#' 
#' @returns 
#' The \link[methods]{show} method for \linkS4class{moment} object 
#' does not have a returned value.
#' 
#' 
#' @export
setMethod(f = show, signature = signature(object = 'moment'), definition = function(object) {
  cat(sprintf('Distribution: %s\n', object@distname))
  cat(sprintf('Mean: %s\n', paste0(sprintf(fmt = '%.3f', object@mean), collapse = ', ')))
  cat(sprintf('Standard Deviations: %s\n', paste0(sprintf(fmt = '%.3f', object@sd), collapse = ', ')))
  cat(sprintf('Skewness: %s\n', paste0(sprintf(fmt = '%.3f', object@skewness), collapse = ', ')))
  cat(sprintf('Kurtosis: %s\n', paste0(sprintf(fmt = '%.3f', object@kurtosis), collapse = ', ')))
})




