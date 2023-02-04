



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
#' @return
#' \link{as.fmx} returns an \linkS4class{fmx} object.
#' 
#' @seealso \link{as.fmx.fitdist} \link{as.fmx.mixEM}
#' 
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
#' @return 
#' \link{as.fmx.fitdist} returns an \linkS4class{fmx} object.
#' 
#' @export
as.fmx.fitdist <- function(x, data = x[['data']], ...) {
  if (!length(data)) stop('Rerun ?fitdistrplus::fitdist with `keepdata = TRUE')
  new(Class = 'fmx', 
      pars = matrix(x[['estimate']], nrow = 1L), distname = x[['distname']],
      data = data, epdf = approxdens(data), 
      vcov = x[['vcov']])
}




#' @title Convert \code{mixEM} Objects to \linkS4class{fmx} Objects
#' 
#' @description ..
#' 
#' @param x \code{mixEM} object
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param ... ..
#' 
#' @note 
#' \link[mixtools]{plot.mixEM} not plot \link[mixtools]{gammamixEM} returns, as of 2022-09-19.
#' 
#' @return 
#' \link{as.fmx.mixEM} returns an \linkS4class{fmx} object.
#' 
#' @examples 
#' library(mixtools)
#' (x = as.fmx(normalmixEM(faithful$waiting, k = 2)))
#' 
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









#' @title Convert Objects from \pkg{mixsmsn} to \linkS4class{fmx} Objects
#' 
#' @description ..
#' 
#' @param x \code{Normal}, \code{Skew.normal}, \code{t}, \code{Skew.t} object, returned from \link[mixsmsn]{smsn.mix} with option \code{family = 'Skew.normal'}
#' 
#' @param data \link[base]{numeric} \link[base]{vector}
#' 
#' @param engine \link[base]{character} scalar for \code{Skew.normal} and \code{Skew.t} method dispatch, 
#' either \CRANpkg{sn} 
#' (default, using \link[sn]{dsn}, \link[sn]{psn}, \link[sn]{dst} or \link[sn]{pst}),
#' or \CRANpkg{mixsmsn} (using my re-written \link{dSN} and \link{dST}, but I do not have function \code{pSN} or \code{pST}).
#'  
#' @param ... ..
#' 
#' @seealso \link[mixsmsn]{smsn.mix}
#' 
#' @note
#' \link[mixsmsn]{smsn.mix} does not offer a parameter to keep the input data, as of 2021-10-06.
#' 
#' @return 
#' 
#' \link{as.fmx.Skew.normal}, \code{as.fmx.Normal} and \link{as.fmx.Skew.t} return an \linkS4class{fmx} object.
#' 
#' \link{as.fmx.t} has not been completed yet
#' 
#' @examples 
#' 
#' if (FALSE) {
#' # all examples work
#' # disabled as they are slow
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
#' class(m <- smsn.mix(x, nu = 3, g = 3, family = 'Skew.normal', calc.im = FALSE))
#' mix.hist(y = x, model = m)
#' 
#' m2 = as.fmx(m, data = x, engine = 'mixsmsn')
#' autoplot(m2)
#' (l2 = logLik(m2))
#' stopifnot(identical(AIC(m2), m$aic), identical(BIC(m2), m$bic))
#' tryCatch(fmx_diagnosis(m2, type = 'Kolmogorov'), error = identity)
#' tryCatch(fmx_diagnosis(m2, type = 'CramerVonMises'), error = identity)
#' 
#' m3 = as.fmx(m, data = x, engine = 'sn')
#' autoplot(m3)
#' (l3 = logLik(m3))
#' stopifnot(identical(AIC(m3), m$aic), identical(BIC(m3), m$bic))
#' autoplot(m3, type = 'distribution')
#' fmx_diagnosis(m3, type = 'Kolmogorov')
#' fmx_diagnosis(m3, type = 'CramerVonMises')
#' 
#' fmx_diagnosis(m2, type = 'KullbackLeibler') - fmx_diagnosis(m3, type = 'KullbackLeibler')
#' 
#' identical(l2, l3)
#' range(attr(l2, 'logl') - attr(l3, 'logl'))
#' 
#' if (FALSE) {
#' (x0 = x[729]) # only wrong
#' #range(qfmx(pfmx(x[-729], dist = m3), dist = m3) - x[-729])
#' #range(qfmx(pfmx(x[-729], dist = m3, lower.tail = F), dist = m3, lower.tail = F) - x[-729])
#' qfmx(pfmx(x0, dist = m3), dist = m3)
#' qfmx(pfmx(x0, dist = m3, lower.tail = FALSE), dist = m3, lower.tail = FALSE)
#' qfmx_interval(m3) # using interval end.   I have no plan to fix this error.
#' }
#' 
#' # Normal
#' class(m <- smsn.mix(x, nu = 3, g = 3, family = 'Normal', calc.im = FALSE))
#' mix.hist(y = x, model = m)
#' autoplot(as.fmx(m, data = x))
#' 
#' 
#' # skew t
#' set.seed(131534); x = rmix(n = 1e3L, p = c(.5, .2, .3), family = 'Skew.t', 
#'   arg = list(unname(arg1), unname(arg2), unname(arg3)))
#'
#' g = 3
#' class(m <- smsn.mix(x, nu = 3, g = g, family = 'Skew.t', calc.im = FALSE))
#' mix.hist(y = x, model = m)
#' 
#' m2 = as.fmx(m, data = x, engine = 'mixsmsn')
#' autoplot(m2)
#' (l2 = logLik(m2))
#' attr(l2, 'df') = g*3 + 1 + (g-1) # ?mixsmsn::smsn.mix gives same `nu` for all components
#' stopifnot(all.equal.numeric(AIC(l2), m$aic), all.equal.numeric(BIC(l2), m$bic))
#' tryCatch(fmx_diagnosis(m2, type = 'Kolmogorov'), error = identity)
#' tryCatch(fmx_diagnosis(m2, type = 'CramerVonMises'), error = identity)
#' 
#' m3 = as.fmx(m, data = x, engine = 'sn')
#' autoplot(m3)
#' (l3 = logLik(m3))
#' attr(l3, 'df') = g*3 + 1 + (g-1) # same reason
#' stopifnot(all.equal.numeric(AIC(l3), m$aic), all.equal.numeric(BIC(l3), m$bic))
#' autoplot(m3, type = 'distribution')
#' fmx_diagnosis(m3, type = 'Kolmogorov')
#' fmx_diagnosis(m3, type = 'CramerVonMises')
#' 
#' fmx_diagnosis(m2, type = 'KullbackLeibler') - fmx_diagnosis(m3, type = 'KullbackLeibler')
#' 
#' identical(l2, l3)
#' range(attr(l2, 'logl') - attr(l3, 'logl'))
#' 
#' # t
#' class(m <- smsn.mix(x, nu = 3, g = 3, family = 't', calc.im = FALSE))
#' mix.hist(y = x, model = m)
#' # autoplot(as.fmx(m, data = x)) # not ready yet!!
#' }
#' 
#' 
#' 
#' @name as_fmx_mixsmsn
#' @method as.fmx Skew.normal
#' @export
as.fmx.Skew.normal <- function(x, data, engine = c('sn', 'mixsmsn'), ...) {
  if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
  x <- sort.Skew.normal(x, decreasing = FALSE)
  engine <- match.arg(engine)
  
  new(Class = 'fmx', pars = switch(engine, sn = cbind(
    xi = x[['mu']],
    omega = sqrt(x[['sigma2']]),
    alpha = x[['shape']]
  ), mixsmsn = cbind(
    mean = x[['mu']],
    sd = sqrt(x[['sigma2']]),
    shape = x[['shape']]
  )), 
  w = x[['pii']],
  distname = switch(engine, sn = 'sn', mixsmsn = 'SN'), 
  data = data, epdf = approxdens(data))
}




#' @rdname as_fmx_mixsmsn
#' @method as.fmx Normal
#' @export
as.fmx.Normal <- function(x, data, ...) {
  if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
  x <- sort.Normal(x, decreasing = FALSE)
  
  new(Class = 'fmx', pars = cbind(
    mean = x[['mu']],
    sd = sqrt(x[['sigma2']])
  ), 
  w = x[['pii']],
  distname = 'norm', 
  data = data, epdf = approxdens(data))
}





#' @rdname as_fmx_mixsmsn
#' @method as.fmx Skew.t
#' @export
as.fmx.Skew.t <- function(x, data, engine = c('sn', 'mixsmsn'), ...) {
  if (!length(data) || !is.numeric(data) || anyNA(data)) stop('illegal `data`')
  x <- sort.Skew.t(x, decreasing = FALSE)
  engine <- match.arg(engine)
  
  new(Class = 'fmx', pars = switch(engine, sn = cbind(
    xi = x[['mu']],
    omega = sqrt(x[['sigma2']]),
    alpha = x[['shape']],
    nu = x[['nu']]
  ), mixsmsn = cbind(
    mean = x[['mu']],
    sd = sqrt(x[['sigma2']]),
    shape = x[['shape']],
    nu = x[['nu']]
  )), 
  w = x[['pii']],
  distname = switch(engine, sn = 'st', mixsmsn = 'ST'), 
  data = data, epdf = approxdens(data))
}





#' @rdname as_fmx_mixsmsn
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
