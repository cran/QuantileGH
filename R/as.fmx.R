



#' @title Turn Various Objects to \linkS4class{fmx}
#' 
#' @description Turn various objects that are created in other packages 
#' to \linkS4class{fmx} class
#' 
#' @param x an R object
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @details 
#' In order to take advantage of all methods for \linkS4class{fmx} objects
#' 
#' @returns
#' S3 generic function [as.fmx()] returns an \linkS4class{fmx} object.
#' 
#' @seealso [as.fmx.fitdist()] [as.fmx.mixEM()]
#' @export
as.fmx <- function(x, data, ...) UseMethod('as.fmx')


#' @export
as.fmx.fmx <- function(x, data, ...) x


#' @title Convert \link[fitdistrplus]{fitdist} Objects to \linkS4class{fmx} Objects
#' 
#' @description ..
#' 
#' @param x \link[fitdistrplus]{fitdist} object
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @returns 
#' Function [as.fmx.fitdist()] returns an \linkS4class{fmx} object.
#' 
#' @export as.fmx.fitdist
#' @export
as.fmx.fitdist <- function(x, data = x[['data']], ...) {
  if (!length(data)) stop('Rerun ?fitdistrplus::fitdist with `keepdata = TRUE')
  new(Class = 'fmx', 
      pars = matrix(x[['estimate']], nrow = 1L), distname = x[['distname']],
      data = data, epdf = approxdens(data), 
      vcov = x[['vcov']])
}




#' @title Convert `mixEM` Objects to \linkS4class{fmx} Objects
#' 
#' @description ..
#' 
#' @param x `mixEM` object
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @note 
#' \link[mixtools]{plot.mixEM} not plot \link[mixtools]{gammamixEM} returns, as of 2022-09-19.
#' 
#' @returns 
#' Function [as.fmx.mixEM()] returns an \linkS4class{fmx} object.
#' 
#' @examples 
#' library(mixtools)
#' (x = as.fmx(normalmixEM(faithful$waiting, k = 2)))
#' 
#' @export as.fmx.mixEM
#' @export
as.fmx.mixEM <- function(x, data = x[['x']], ...) {
  if (!length(data)) stop('wont happen')
  x <- sort.mixEM(x, decreasing = FALSE)
  
  switch(x[['ft']], normalmixEM = {
    pars <- cbind(mean = x[['mu']], sd = x[['sigma']])
    distname <- 'norm'
  }, gammamixEM = {
    pars <- t.default(x[['gamma.pars']])
    colnames(pars) <- c('shape', 'scale') # names of parameters of ?stats::dgamma
    # read \link[mixtools]{gammamixEM} carefully: 'beta' is actually `scale`
    distname <- 'gamma'
  }, stop(x[['ft']], ' not supported yet'))
  new(Class = 'fmx', 
      pars = pars, w = x[['lambda']], distname = distname,
      data = data, epdf = approxdens(data))
}









#' @title Convert `Skew.normal` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description ..
#' 
#' @param x `'Skew.normal'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 'Skew.normal'`
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.Skew.normal()] returns an \linkS4class{fmx} object.
#' 
#' @examples 
#' 
#' library(mixsmsn)
#' # ?smsn.mix
#' arg1 = c(mu = 5, sigma2 = 9, lambda = 5, nu = 5)
#' arg2 = c(mu = 20, sigma2 = 16, lambda = -3, nu = 5)
#' arg3 = c(mu = 35, sigma2 = 9, lambda = -6, nu = 5)
#' set.seed(120); x = rmix(n = 1e3L, p=c(.5, .2, .3), family = 'Skew.t', 
#'   arg = list(unname(arg1), unname(arg2), unname(arg3)))
#'
#' # Skew Normal
#' class(m1 <- smsn.mix(x, nu = 3, g = 3, family = 'Skew.normal', calc.im = FALSE))
#' mix.hist(y = x, model = m1)
#' m1a = as.fmx(m1, data = x)
#' autoplot(m1a)
#' (l1a = logLik(m1a))
#' stopifnot(identical(AIC(m1a), m1$aic), identical(BIC(m1a), m1$bic))
#' autoplot(m1a, type = 'distribution')
#' if (FALSE) {
#' range(qfmx(pfmx(x, dist = m1a), dist = m1a) - x) # need to think why
#' # may have to do with ?qfmx_interval
#' }
#' 
#' 
#' @method as.fmx Skew.normal
#' @export as.fmx.Skew.normal
#' @export
as.fmx.Skew.normal <- function(x, data, ...) {
  x <- sort.Skew.normal(x, decreasing = FALSE)
  
  ret <- new(Class = 'fmx', pars = cbind(
    xi = x[['mu']],
    omega = sqrt(x[['sigma2']]),
    alpha = x[['shape']]
  ), 
  w = x[['pii']],
  distname = 'sn')
  
  if (!missing(data)) {
    if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
    ret@data <- data
    ret@epdf <- approxdens(data)
    ret@Kolmogorov <- Kolmogorov_fmx(ret)
    ret@CramerVonMises <- CramerVonMises_fmx(ret)
    ret@KullbackLeibler <- KullbackLeibler_fmx(ret)
  }
  
  return(ret)
}





#' @title Convert `Normal` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description ..
#' 
#' @param x `'Normal'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 'Normal'`
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.Normal()] returns an \linkS4class{fmx} object.
#' 
#' @examples 
#' 
#' library(mixsmsn)
#' # ?smsn.mix
#' arg1 = c(mu = 5, sigma2 = 9, lambda = 5, nu = 5)
#' arg2 = c(mu = 20, sigma2 = 16, lambda = -3, nu = 5)
#' arg3 = c(mu = 35, sigma2 = 9, lambda = -6, nu = 5)
#' set.seed(120); x = rmix(n = 1e3L, p=c(.5, .2, .3), family = 'Skew.t', 
#'   arg = list(unname(arg1), unname(arg2), unname(arg3)))
#'
#' # Normal
#' class(m2 <- smsn.mix(x, nu = 3, g = 3, family = 'Normal', calc.im = FALSE))
#' mix.hist(y = x, model = m2)
#' m2a = as.fmx(m2, data = x)
#' autoplot(m2a)
#' 
#' @method as.fmx Normal
#' @export as.fmx.Normal
#' @export
as.fmx.Normal <- function(x, data, ...) {
  x <- sort.Normal(x, decreasing = FALSE)
  
  ret <- new(Class = 'fmx', pars = cbind(
    mean = x[['mu']],
    sd = sqrt(x[['sigma2']])
  ), 
  w = x[['pii']],
  distname = 'norm')
  
  if (!missing(data)) {
    if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
    ret@data <- data
    ret@epdf <- approxdens(data)
    ret@Kolmogorov <- Kolmogorov_fmx(ret)
    ret@CramerVonMises <- CramerVonMises_fmx(ret)
    ret@KullbackLeibler <- KullbackLeibler_fmx(ret)
  }
  
  return(ret)
  
}






#' @title Convert `Skew.t` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description ..
#' 
#' @param x `'Skew.t'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 'Skew.t'` 
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.Skew.t()] returns an \linkS4class{fmx} object.
#' 
#' @examples 
#' \donttest{
#' # mixsmsn::smsn.mix with option `family = 'Skew.t'` is slow
#' 
#' library(mixsmsn)
#' # ?smsn.mix
#' arg1 = c(mu = 5, sigma2 = 9, lambda = 5, nu = 5)
#' arg2 = c(mu = 20, sigma2 = 16, lambda = -3, nu = 5)
#' arg3 = c(mu = 35, sigma2 = 9, lambda = -6, nu = 5)
#' set.seed(120); x = rmix(n = 1e3L, p=c(.5, .2, .3), family = 'Skew.t', 
#'   arg = list(unname(arg1), unname(arg2), unname(arg3)))
#'
#' # Skew t
#' class(m3 <- smsn.mix(x, nu = 3, g = 3, family = 'Skew.t', calc.im = FALSE))
#' mix.hist(y = x, model = m3)
#' m3a = as.fmx(m3, data = x)
#' autoplot(m3a)
#' (l3a = logLik(m3a))
#' stopifnot(all.equal.numeric(AIC(l3a), m3$aic), all.equal.numeric(BIC(l3a), m3$bic))
#' autoplot(m3a, type = 'distribution')
#' }
#' 
#' @method as.fmx Skew.t
#' @export as.fmx.Skew.t
#' @export
as.fmx.Skew.t <- function(x, data, ...) {
  x <- sort.Skew.t(x, decreasing = FALSE)
  
  K <- length(x[['mu']]) # number of components
  if (length(x[['nu']]) != 1L) stop('\\pkg{mixsmsn} update to enable multiple `nu`? Modify ?npar.fmx') 
  
  ret <- new(Class = 'fmx', pars = cbind(
    xi = x[['mu']],
    omega = sqrt(x[['sigma2']]),
    alpha = x[['shape']],
    nu = x[['nu']]
  ), 
  w = x[['pii']],
  distname = 'st')
  
  if (!missing(data)) {
    if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
    ret@data <- data
    ret@epdf <- approxdens(data)
    ret@Kolmogorov <- Kolmogorov_fmx(ret)
    ret@CramerVonMises <- CramerVonMises_fmx(ret)
    ret@KullbackLeibler <- KullbackLeibler_fmx(ret)
  }
  
  return(ret)
}






#' @title Convert `Normal` fit from \CRANpkg{mixsmsn} to \linkS4class{fmx}
#' 
#' @description ..
#' 
#' @param x `'t'` object, 
#' returned from \link[mixsmsn]{smsn.mix} with parameter 
#' `family = 't'`
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @returns 
#' Function [as.fmx.t()] has not been completed yet
#' 
#' @examples 
#' 
#' library(mixsmsn)
#' # ?smsn.mix
#' arg1 = c(mu = 5, sigma2 = 9, lambda = 5, nu = 5)
#' arg2 = c(mu = 20, sigma2 = 16, lambda = -3, nu = 5)
#' arg3 = c(mu = 35, sigma2 = 9, lambda = -6, nu = 5)
#' set.seed(120); x = rmix(n = 1e3L, p=c(.5, .2, .3), family = 'Skew.t', 
#'   arg = list(unname(arg1), unname(arg2), unname(arg3)))
#'
#' # t
#' class(m4 <- smsn.mix(x, nu = 3, g = 3, family = 't', calc.im = FALSE))
#' mix.hist(y = x, model = m4)
#' # autoplot(as.fmx(m4, data = x)) # not ready yet!!
#' 
#' 
#' @method as.fmx t
#' @export
as.fmx.t <- function(x, data, ...) {
  if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
  x <- sort.t(x, decreasing = FALSE)
  stop('need to programe scale_and_shift_t')
  #new(Class = 'fmx', pars = cbind(
  #  mean = x[['mu']],
  #  sd = sqrt(x[['sigma2']])
  #), 
  #w = x[['pii']],
  #distname = 'norm', 
  #data = data, epdf = approxdens(data))
}
