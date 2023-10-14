


#' @title Show \linkS4class{fmx} Object
#' 
#' @description
#' Print the parameters of an \linkS4class{fmx} object and plot its density curves.
#' 
#' @param object an \linkS4class{fmx} object
#' 
#' @returns 
#' The \link[methods]{show} method for \linkS4class{fmx} object 
#' does not have a returned value.
#' 
#' @export
setMethod(f = show, signature = signature(object = 'fmx'), definition = function(object) {
  print.fmx(object)
  if (length(test_table <- attr(object, which = 'steps', exact = TRUE))) {
    direction <- attr(object, which = 'direction', exact = TRUE)
    cat('Stepwise Parameter', switch(direction, drop1 = 'Elimination', add1 = 'Growth'), '(Fixed # of Components)\n')
    print.data.frame(test_table)
    cat('\n')
  }
  print(autoplot.fmx(object))
  return(invisible())
})




#' @title S3 \link[base]{print} of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param x an \linkS4class{fmx} object
#' 
#' @param ... additional parameters, not currently in use
#' 
#' @returns 
#' [print.fmx] returns the input \linkS4class{fmx} object invisibly.
#' 
#' @seealso 
#' \link[base]{print}
#' 
#' @export
print.fmx <- function(x, ...) {
  pars <- x@pars
  K <- dim(pars)[1L]
  pars[] <- sprintf(fmt = '%.2f', pars)
  dimnames(pars) <- list(paste0(seq_len(K), '-comp.'), distArgs(x@distname))
  if (length(id_constr <- fmx_constraint(x))) pars[id_constr] <- '.'
  obj <- if (K == 1L) pars else cbind(pars, w = sprintf(fmt = '%.1f%%', x@w*1e2))
  heading <- paste0(K, '-Component Mixture of ', sQuote(x@distname), ' Distribution')
  
  if (length(x@data)) {
    # theoretically does not make sense to talk about confidence interval, without a sample
    ci <- confint.fmx(x, level = .95, internal = FALSE)
    id_constr <- fmx_constraint(x)
    if (length(ci) && !anyNA(ci)) {
      ci0 <- sprintf(fmt = '(%.2f~%.2f)', ci[,1L], ci[,2L])
      if (length(id_constr)) {
        obj[id_constr] <- '.'
        obj[-id_constr] <- paste(obj[-id_constr], ci0)
      } else obj[] <- paste(obj, ci0)
      heading <- paste0(heading,  ' (w. 95% Confidence Intervals)')
    } #else heading <- paste0('Malformed ', heading)
  }
  
  cat('\n ', heading, '\n\n', sep = '')
  print.default(obj, quote = FALSE)
  if (length(id_constr)) cat('\nwhere ', sQuote('.'), ' denotes an enforced constraint\n', sep = '')
  cat('\n')
  
  return(invisible(x))
}







#' @title Subset of Components in \linkS4class{fmx} Object
#' 
#' @description 
#' 
#' Taking subset of components in \linkS4class{fmx} object
#' 
#' @param x \linkS4class{fmx} object
#' 
#' @param i \link[base]{integer} or \link[base]{logical} \link[base]{vector}, 
#' the row indices of *components* to be chosen, see \link[base]{[}
#' 
#' @details 
#' 
#' Note that using definitions as S3 method dispatch \code{`[.fmx`} won't work 
#' for \linkS4class{fmx} objects.
#' 
#' @returns 
#' 
#' An \linkS4class{fmx} object consisting of a subset of components.
#' information about the observations (e.g. slots `@@data` and `@@data.name`),
# as well as other estimation related slots (e.g., `@@init`) 
#' will be lost.
#' 
#' @examples 
#' 
#' (d = fmx('norm', mean = c(1, 4, 7), w = c(1, 1, 1)))
#' d[1:2]
#' 
#' @export
setMethod(`[`, signature(x = 'fmx', i = 'ANY'), definition = function(x, i) {
  if (length(x@data)) message('Subset the estimates and drop `@data` etc.')
  pars <- x@pars[i, , drop = FALSE]
  w <- x@w[i]
  w <- unname(w / sum(w)) # adjust mixing proportions
  o <- order(pars[, 1L])
  new(Class = 'fmx', pars = pars[o, , drop = FALSE], w = w[o], distname = x@distname)
})










#' @title Number of Observations in \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' [nobs.fmx] returns the sample size of
#' the observations used in [QLMDe] estimation, or `integer(0)` for distribution-only
#' \linkS4class{fmx} object 
#' 
#' @returns 
#' 
#' [nobs.fmx] returns an \link[base]{integer} scalar.
#' 
#' @seealso 
#' \link[stats]{nobs} 
#' 
#' @export nobs.fmx
#' @export
nobs.fmx <- function(object, ...) {
  if (length(n <- object@data)) return(n)
  return(integer())
}








#' @title Confidence Interval of \linkS4class{fmx} Object
#' 
#' @description ...
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param level confidence level, default \eqn{95\%}.
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' [confint.fmx] returns the Wald-type confidence intervals based on the user-friendly parameters (`parm = 'user'`),
#'  or the internal/unconstrained parameters (`parm = 'internal'`).
#' When the distribution has constraints on one or more parameters, 
#' \link{confint.fmx} does not return the confident intervals of for the constrained parameters.
#'  
#' @returns 
#' [confint.fmx] returns a \link[base]{matrix}
#' 
#' @importFrom stats confint
#' @export confint.fmx
#' @export
confint.fmx <- function(object, ..., level = .95) {
  # essentially ?stats::confint.default
  cf <- coef.fmx(object, ...)
  if (!length(vv <- vcov.fmx(object, ...))) return(invisible())
  ses <- sqrt(diag(vv))
  p1 <- (1 - level) / 2
  p <- c(p1, 1 - p1)
  ret <- cf + ses %*% t.default(qnorm(p))
  dimnames(ret) <- list(names(cf), sprintf('%.1f%%', 1e2*p))
  return(ret)
}



#' @title Variance-Covariance of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param internal \link[base]{logical} scalar, either for the user-friendly parameters (`FALSE`, default)
#' (e.g., `mean,sd` for normal mixture, and `A,B,g,h` for Tukey's \eqn{g}-and-\eqn{h} mixture), or
#' for the internal/unconstrained parameters (`TRUE`).
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' [vcov.fmx] returns 
#' the approximate asymptotic variance-covariance \link[base]{matrix} of the user-friendly parameters via delta-method (`parm = 'user'`), 
#' or the asymptotic variance-covariance matrix of the internal/unconstrained parameters (`parm = 'internal'`). 
#' When the distribution has constraints on one or more parameters, 
#' [vcov.fmx] does not return the variance/covariance involving the constrained parameters.
#' 
#' @returns 
#' 
#' [vcov.fmx] returns a \link[base]{matrix}.
#' 
#' @importFrom stats vcov
#' @export vcov.fmx
#' @export
vcov.fmx <- function(object, internal = FALSE, ...) {
  
  if (!length(object@data)) return(invisible())
  
  if (!internal && length(vv <- object@vcov)) return(vv) # 'fitdist' objects
  
  int_vv <- object@vcov_internal
  if (internal) return(int_vv)
  
  if (!length(int_vv)) return(int_vv) # wont be able to computer user-vcov if internal-vcov is wrong
  
  distname <- object@distname
  pars <- object@pars
  K <- dim(pars)[1L]
  int_nm <- dimnames(int_vv)[[1L]]
  
  int_p <- fmx2dbl(object) # internal parameters
  anm <- distArgs(distname)
  n_anm <- length(anm)
  user_nm <- c(t.default(outer(c(anm, if (K > 1L) 'w'), 1:K, FUN = paste0)))
  jacob <- array(0, dim = c(length(int_p), length(user_nm)), dimnames = list(names(int_p), user_nm))
  # location parameters A_1 -> A_1
  jacob[1L, 1:K] <- 1
  if (K > 1L) {
    # location parameters A_2, .., A_k -> d_2, .., d_k
    for (k in 2:K) jacob[k, k:K] <- exp(int_p[k]) 
    # mixture parameters w_1, .., w_k -> pi_2, .., pi_k
    id_pi <- (n_anm*K+1L):((n_anm+1L)*K-1L)
    e_pi <- exp(int_p[id_pi])
    sum_pi <- sum(1 + e_pi)^2 # 1 = exp(pi_1) = exp(0)
    jacob[id_pi, n_anm*K+1L] <- - e_pi / sum_pi
    jacob[id_pi, id_pi+1L] <- - tcrossprod(e_pi) / sum_pi
  }
  
  switch(distname, norm = {
    id_exp <- (K+1L):(2*K) # 'sd'
    id_identity <- id_constr <- NULL
  }, GH = {
    id_exp <- c((K+1L):(2*K), (3*K+1L):(4*K)) # 'B' and 'h'
    id_identity <- (2*K+1L):(3*K) # 'g'
    id_constr <- fmx_constraint(object)
  })
  
  if (length(id_exp)) jacob[cbind(id_exp, id_exp)] <- exp(int_p[id_exp])
  if (length(id_identity)) jacob[cbind(id_identity, id_identity)] <- 1
  jacob_free <- if (length(id_constr)) jacob[int_nm, -id_constr] else jacob
  
  return(t.default(jacob_free) %*% int_vv %*% jacob_free)
  
}





# ?stats::residuals 
# @export
# residuals.fmx <- function(object, ...) stop('useful?')





#' @title Parameter Estimates of \linkS4class{fmx} object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param internal \link[base]{logical} scalar, either for the user-friendly parameters (`FALSE`, default)
#' (e.g., `mean,sd` for normal mixture, and `A,B,g,h` for Tukey's \eqn{g}-and-\eqn{h} mixture), or
#' for the internal/unconstrained parameters (`TRUE`).
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' Function [coef.fmx()] returns the estimates of the user-friendly parameters (`parm = 'user'`), 
#' or the internal/unconstrained parameters (\code{parm = 'internal'}).
#' When the distribution has constraints on one or more parameters, 
#' [coef.fmx()] does not return the estimates (which is constant \code{0}) of the constrained parameters.
#' 
#' @returns 
#' 
#' Function [coef.fmx()] returns a \link[base]{numeric} \link[base]{vector}.
#' 
#' @importFrom stats coef
#' @export coef.fmx
#' @export
coef.fmx <- function(object, internal = FALSE, ...) {
  anm <- distArgs(object@distname)
  K <- dim(object@pars)[1L]
  cf0 <- if (internal) fmx2dbl(object) else if (K == 1L) {
    setNames(c(object@pars), nm = c(t.default(outer(anm, 1:K, FUN = paste0))))
  } else {
    setNames(c(object@pars, object@w), nm = c(t.default(outer(c(anm, 'w'), 1:K, FUN = paste0))))
  }
  if (!length(id_constr <- fmx_constraint(object))) return(cf0)
  return(cf0[-id_constr])
}






#' @title Log-Likelihood of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param data \link[base]{double} \link[base]{vector}, actual observations
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' [logLik.fmx] returns a \link[stats]{logLik} object indicating the log-likelihood.
#' An additional attribute `attr(,'logl')` indicates the point-wise log-likelihood, 
#' to be use in Vuong's closeness test.
#' 
#' @returns 
#' 
#' [logLik.fmx] returns a \link[stats]{logLik} object with 
#' an additional attribute `attr(,'logl')`.
#' 
#' @importFrom stats logLik
#' @export logLik.fmx
#' @export
logLik.fmx <- function(object, data = object@data, ...) {
  
  # for developer to batch-calculate AIC/BIC quickly
  #if (length(objF <- attr(object, which = 'objF', exact = TRUE))) {
  #  if (inherits(objF[[1L]], what = 'logLik')) return(objF[[1L]])
  #}
  # ?step_fmx no longer uses logLik
  
  if (!length(data)) return(invisible())
  
  logd <- dfmx(x = data, dist = object, log = TRUE, ...)
  if (!all(is.finite(logd))) {
    #print(logd)
    #object <<- object
    #stop('malformed fit (?.dGH has been well debug-ged)')
    # very likely to be `B = 0`
    # do not stop.  settle with -Inf log-likelihood
  }
  ret <- sum(logd)
  attr(ret, which = 'logl') <- logd # additional attributes; needed in Vuong's test
  attr(ret, which = 'nobs') <- length(data)
  attr(ret, which = 'df') <- npar.fmx(object)
  class(ret) <- 'logLik'
  return(ret)
}






#' @title Create \link[ggplot2]{layer} for Continuous \linkS4class{fmx} Objects
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param type \link[base]{character} scalar.  
#' Option `'density'` (default) plots the probability density for \linkS4class{fmx} input
#' (and the histogram if argument `data` is available).
#' Option `'distribution'` plots the cumulative probability distribution for \linkS4class{fmx} input 
#' (and the empirical cumulative distribution if argument `data` is available).
#' 
#' @param data (optional) \link[base]{numeric} \link[base]{vector} of the observations.
#' Default is the slot `object@@data`.
#' 
#' @param epdf (optional) empirical probability density \link[base]{function} returned by \link[stats]{approxfun}.
#' Default is the slot `object@@epdf`
#' 
#' @param hist.fill color of the body of histogram, default `'grey95'`
#' 
#' @param n \link[base]{integer}, see \link[ggplot2]{stat_function}
#' 
#' @param curve.col color of the density curve of the fitted finite mixture distribution.
#' Default `'black'`
#' 
#' @param xlim \link[base]{numeric} length-two \link[base]{vector}, horizontal range
#' 
# @param init \link[base]{logical} scalar, whether to plot the initial estimates used in \link{QLMDe}, default `FALSE`.
#' 
#' @param probs \link[base]{numeric} \link[base]{vector}, 
#' the percentages (to be) used in \link{QLMDe}, can be plotted as vertical lines.
#' Use `probs = NULL` to suppress the printing of these lines.
#' 
#' @param ... potential parameters of \link[ggplot2]{stat_function}
#' 
#' @returns 
#' 
#' Function [autolayer_fmx_continuous()] returns a \link[base]{list} of \link[ggplot2]{layer}s.
#' 
#' @seealso 
#' \link[ggplot2]{autolayer} 
#' 
#' @importFrom ggplot2 geom_histogram stat_function stat_ecdf geom_vline geom_text
#' @importFrom scales percent
#' @export
autolayer_fmx_continuous <- function(
    object, 
    type = c('density', 'distribution'), 
    data = object@data, epdf = object@epdf,
    probs = object@probs,
    # init = attr(object, which = 'init', exact = TRUE), 
    xlim = if (!length(data)) qfmx(p = c(.01, .99), dist = object) else range.default(data),
    hist.fill = 'grey95', 
    curve.col = 1, # black curve, default
    n = 1001L,
    ...
) {
  
  type <- match.arg(type)
  fun <- switch(type, density = dfmx, distribution = pfmx)
  
  data_lyr <- if (length(data)) {
    switch(type, density = list(
      geom_histogram(mapping = aes(x = data, y = after_stat(density)), colour = 'white', fill = hist.fill, bins = 30L),
      if (is.function(epdf)) stat_function(fun = epdf, n = n, colour = 'grey70', linetype = 2L, ...)
    ), distribution = list(
      stat_ecdf(mapping = aes(x = data), geom = 'step', pad = FALSE, colour = 'grey70', linetype = 2L)
    ))
  } # else NULL
  
  prob_lyr <- if (length(probs)) {
    v <- if (length(data)) quantile(x = data, probs = probs) else qfmx(p = probs, dist = object)
    xlim <- range.default(xlim, v)
    vnm <- sprintf('%.1f%%', 1e2*probs)
    list(
      geom_vline(mapping = aes(xintercept = v, colour = vnm), size = .1, show.legend = FALSE),
      geom_text(mapping = aes(x = v, y = 0, colour = vnm, label = vnm, hjust = -.5), size = 3, angle = 90, show.legend = FALSE)
    )
  } # else NULL
  
  ret <- list(
    data_lyr, prob_lyr,
    stat_function(fun = fun, args = list(dist = object), xlim = xlim, n = n, colour = curve.col, ...),
    switch(type, distribution = scale_y_continuous(labels = percent))
    # if (inherits(init, what = 'fmx')) stat_function(fun = fun, args = list(dist = init), n = n, linetype = 2L, colour = curve.col, ...),
  )
  return(ret[lengths(ret) > 0L])
  
}








#' @title Create \link[ggplot2]{layer} for Discrete \linkS4class{fmx} Objects
#' 
#' @description ..
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param type \link[base]{character} scalar.  
#' Option `'density'` (default) plots the probability density for \linkS4class{fmx} input
#' (and the histogram if argument `data` is available).
#' Option `'distribution'` plots the cumulative probability distribution for \linkS4class{fmx} input 
#' (and the cumulative histogram if argument `data` is available).
#' 
#' @param data (optional) \link[base]{numeric} (actually \link[base]{integer}) \link[base]{vector} of the observations.
#' Default is the slot `object@@data`.
#' 
#' @param bins \link[base]{integer} scalar
#' 
#' @param xlim \link[base]{numeric} length-two \link[base]{vector}, horizontal range
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' 
#' [autolayer_fmx_discrete] returns a \link[base]{list} of \link[ggplot2]{layer}s.
#' 
#' @seealso 
#' \link[ggplot2]{autolayer}
#' 
#' @importFrom ggplot2 geom_histogram stat_ecdf
#' @importFrom scales percent
#' @export
autolayer_fmx_discrete <- function(
    object, 
    type = c('density', 'distribution'), 
    data = object@data,
    xlim = if (length(data)) data else qfmx(p = c(.01, .99), dist = object),
    bins = 60L,
    ...
) {
  
  type <- match.arg(type)
  fun <- switch(type, density = dfmx, distribution = pfmx)
  
  if (anyNA(xlim) || !(nL <- length(xlim))) stop('`xlim` must be integer')
  x0 <- if (nL == 1L) {
    if (xlim <= 0) stop('right end of `xlim` must be positive')
    0L:ceiling(xlim)
  } else {
    if (any(xlim < 0L)) stop('`xlim` must be non-negative')
    floor(min(xlim)):ceiling(max(xlim))
  } 
  
  y <- fun(x0, dist = object)
  
  return(list(
    if (length(data)) {
      switch(type, density = {
        geom_histogram(mapping = aes(x = data, y = after_stat(density)), bins = bins, colour = 'white', alpha = .1, show.legend = FALSE)
      }, distribution = {
        stat_ecdf(mapping = aes(x = data), geom = 'step', pad = FALSE, colour = 'grey70', linetype = 2L, show.legend = FALSE)
      })
    },
    geom_step(mapping = aes(x = x0, y = y), stat = 'identity'),
    scale_y_continuous(labels = percent)
  ))
  
}










#' @title Plot \linkS4class{fmx} Objects using \CRANpkg{ggplot2}
#' 
#' @description 
#' 
#' Plot \linkS4class{fmx} objects using \CRANpkg{ggplot2}.
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param xlab,ylab,title,caption \link[base]{character} scalars, the 
#' horizontal and vertical label, title and caption
#' 
#' @param ... potential parameters of [autolayer_fmx_continuous] and [autolayer_fmx_discrete]
#' 
#' @returns 
#' 
#' [autoplot.fmx] returns a \link[ggplot2]{ggplot} object.
#' 
#' @seealso [autolayer_fmx_continuous] [autolayer_fmx_discrete]
#' \link[ggplot2]{autoplot}
#' 
#' 
#' 
#' @examples 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' curve(dfmx(x, dist = d2), xlim = c(-3, 11))
#' curve(pfmx(x, dist = d2), xlim = c(-3, 11))
#' autoplot(d2)
#' autoplot(d2, type = 'distribution')
#' 
#' @importFrom ggplot2 ggplot labs
#' @importFrom latex2exp TeX
#' @export
autoplot.fmx <- function(
    object, 
    xlab = attr(object, which = 'data.name', exact = TRUE), # S4 slot, 
    ylab = NULL, # paste(object@distname, 'mixture'), 
    title = TeX(getTeX(object)),
    caption = NULL, #if (nv <- length(v)) paste(nv, 'percentiles to match'),
    ...
) {
  ggplot() + 
    (if (object@distname %in% c(distType('continuous'), distType('nonNegContinuous'))) autolayer_fmx_continuous else autolayer_fmx_discrete)(object, ...) +
    labs(x = xlab, y = ylab, title = title, caption = caption)
}





