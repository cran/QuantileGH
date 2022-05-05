




#' @title Specification of \linkS4class{fmx_QLMDe} Class
#' 
#' @description 
#' Quantile least Mahalanobis distance estimates of one-dimensional finite mixture distribution.
#' The \linkS4class{fmx_QLMDe} object contains (i.e., inherits from) the \linkS4class{fmx} object. 
#' 
#' @slot data \link[base]{numeric} vector, the one-dimensional observations
#' 
#' @slot data.name \link[base]{character} scalar, a human-friendly name of observations
#'
#' @slot epdf empirical probability density \link[base]{function} fitted by \link[stats]{approxfun}
#' 
#' @slot quantile_vv variance-covariance \link[base]{matrix} of selected quantiles (based on the selected probabilities stored in slot \code{@@p})
#' 
#' @slot vcov variance-covariance \link[base]{matrix} of the internal (i.e., unconstrained) estimates
#' 
#' @slot init \linkS4class{fmx} object, the initial values to be sent to \link[stats]{optim}
#' 
#' @slot p \link[base]{numeric} vectors of probabilities, where the distance between the empirical and true quantiles are measured
#' 
#' @slot optim a \link[base]{list} returned from \link[stats]{optim}
#' 
#' @export
setClass(Class = 'fmx_QLMDe', contains = 'fmx', slots = c(
  data = 'numeric', 
  data.name = 'character',
  epdf = 'function', 
  quantile_vv = 'matrix',
  vcov = 'matrix',
  init = 'fmx',
  p = 'numeric',
  optim = 'list'
), prototype = prototype(
  data.name = 'observations'
), validity = function(object) {
  if (!length(object@data)) stop('*Never* remove @data slot from \'fmx_QLMDe\' object')
  if (anyNA(object@data)) stop('Observations in \'fmx_QLMDe\' must be free of NA_real_')
  if (!length(object@p)) stop('`p` must be recorded')
})






#' @title Quantile Least Mahalanobis Distance estimates
#' 
#' @description 
#' 
#' The quantile least Mahalanobis distance algorithm estimates the parameters of 
#' single-component or finite mixture distributions   
#' by minimizing the Mahalanobis distance between the vectors of sample and theoretical quantiles.
#' See \link{QLMDp} for the default selection of probabilities at which the sample and theoretical quantiles are compared.
#' 
#' The default initial values are estimated based on trimmed \eqn{k}-means 
#' clustering with re-assignment.
#' 
#' @param x \link[base]{numeric} vector, the one-dimensional observations.
#' 
#' @param data.name \link[base]{character} value, name for the observations for user-friendly print out.
#' 
#' @param distname \link[base]{character} value, name of mixture distribution to be fitted.  Currently supports \code{'norm'} and \code{'GH'}.
#' 
#' @param K \link[base]{integer} scalar, number of components (e.g., must use \code{2L} instead of \code{2}).
#' 
#' @param p \link[base]{numeric} vector, percentiles at where the sample and theoretical quantiles are to be matched.
#' See \code{\link{QLMDp}} for details.
#' 
#' @param init \link[base]{character} scalar for the method of initial values selection, 
#' or an \linkS4class{fmx} object of the initial values. 
#' See \link{QLMDinit} for more details.
#' 
#' @param constraint \link[base]{character} vector, parameters (\eqn{g} and/or \eqn{h} for Tukey's \eqn{g}-&-\eqn{h} mixture) to be set at 0.  
#' See \link{fmx_constraint} for details.
#' 
#' @param tol,maxiter see \link{vuniroot2}
#' 
#' @param ... additional parameters of \code{\link[stats]{optim}}.
#' 
#' @details
#' 
#' Quantile Least Mahalanobis Distance estimator fits a single-component or finite mixture distribution 
#' by minimizing the Mahalanobis distance between
#' the theoretical and observed quantiles,
#' using the empirical quantile variance-covariance matrix \code{\link{quantile_vcov}}.
#' 
#' @return An \linkS4class{fmx_QLMDe} object
#' 
#' @seealso \link[stats]{optim} \link{QLMDinit}
#' 
#' @examples 
#' 
#' \donttest{
#' # Generated from 1-component normal; start with 2-component normal fit
#' set.seed(1623); (y0n <- QLMDe(rnorm(1e3L), distname = 'norm', K = 2L))
#' (y1n <- StepK_fmx(y0n, test = 'logLik', Kmax = 2L)) # one-component
#' vcov(y1n)
#' 
#' # Generated from 2-component normal; start with 1-component normal fit
#' (d1 <- fmx('norm', mean = c(0, 1.5), sd = .5, w = c(.4, .6)))
#' set.seed(100); hist(x1 <- rfmx(n = 1e3L, dist = d1))
#' StepK_fmx(QLMDe(x1, distname = 'norm', K = 1L), test = 'logLik', Kmax = 2L)
#' }
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 1.2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(3123); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' # using user-specified constraint
#' system.time(QLMDe(x2, distname = 'GH', K = 2L, constraint = c('g1', 'h2')))
#' \donttest{
#' # using Step_fmx
#' system.time(y2gh <- QLMDe(x2, distname = 'GH', K = 2L)) # ~2 secs
#' y2gh
#' ks_test(y2gh)
#' CvM_test(y2gh)
#' KL_dist(y2gh)
#' Step_fmx(y2gh, test = 'logLik') # identified true constraint :)
#' }
#' 
#' \donttest{
#' system.time(y1gh <- QLMDe(x2, distname = 'GH', K = 1L))
#' y1gh
#' StepK_fmx(y1gh, test = 'logLik', Kmax = 2L) # correct
#' 
#' #set.seed(1323); x3 <- rGH(n = 1e3L, g = .2, h = .1)
#' #system.time(tmp <- QLMDe(x3, distname = 'GH', K = 2L)) # ~2 secs
#' #StepK_fmx(tmp, Kmax = 2L) # very difficult to drop K ..
#' }
#' 
#' @name QLMDe
#' @export
QLMDe <- function(
  x, distname = c('norm', 'GH'), K, data.name = deparse1(substitute(x)),
  constraint = character(),
  p = QLMDp(x = x),
  init = c('logLik', 'letterValue', 'normix'),
  tol = .Machine$double.eps^.25, maxiter = 1000,
  ...
) {
  
  distname <- match.arg(distname)
  if (length(K) != 1L || !is.numeric(K) || is.na(K) || K <= 0L) stop('number of component must be length-1 positive integer')
  if (!is.integer(K)) stop('number of component must be length-1 positive integer (e.g., use integer `2L` instead of numeric/double `2` for 2-component mixture model)')
  
  if (!is.vector(x, mode = 'double')) stop('x must be double vector')
  if (anyNA(x)) stop('do not allow NA_real_ in observations `x`')
  if (!is.character(data.name) || length(data.name) != 1L || anyNA(data.name) || !nzchar(data.name)) stop('data.name must be length-1 character')
  
  interval <- QLMDe_interval(x = x, ...)
  
  if (anyNA(p)) stop('do not allow NA_real_ in `p`')
  p <- sort.int(unique_allequal(p))
  
  if (missing(init)) {
    init <- switch(match.arg(init), logLik = {
      QLMDinit(x, test = 'logLik', distname = distname, K = K, constraint = constraint)
    }, letterValue = {
      QLMDinit_letterValue(x, distname = distname, K = K, constraint = constraint)
    }, normix = {
      QLMDinit_normix(x, distname = distname, K = K)
    })
  }
  if (!inherits(init, what = 'fmx')) stop('`init` must be \'fmx\' object (e.g., a returned object from ?QLMDinit_letterValue)')
  if ((init@distname != distname) || (dim(init@parM)[1L] != K)) stop('`init` is not a ', distname, ' ', K, '-component fit.')
  
  q_init <- qfmx(p = p, interval = interval, distname = distname, K = K, parM = init@parM, w = init@w)
  if (any(id1 <- is.infinite(q_init))) {
    if (all(id1)) stop('starting values too far away?')
    p <- p[!id1]
  }
  q_obs <- quantile(x, probs = p) # observed quantiles, constant in ?stats::optim
  
  x_kern <- density.default(x) # ?stats::approx inside ?stats::density.default
  x_epdf <- approxfun(x = x_kern$x, y = x_kern$y) # another 'layer' of ?stats::approxfun
  # Tingting is not sure whether ?stats::approx \strong{and} ?stats::approxfun will make `x_epdf(q_obs) = 0` more likely
  d_obs <- x_epdf(q_obs) # observed density evaluated at `q_obs`
  if (anyNA(d_obs)) stop('do not allow NA_real_ empirical density')
  tol <- sqrt(sqrt(.Machine$double.eps))
  if (all(d0 <- (abs(d_obs) < tol))) stop('must have at least one positive density') # `d_obs` should always be positive, but ?base::abs should not hurt
  if (any(d0)) {
    p <- p[!d0] # quantiles, where empirical densities are too close to 0, are excluded from the calculation 
    d_obs <- d_obs[!d0]
    q_obs <- q_obs[!d0]
  }
  
  npar <- K * switch(distname, norm = 2L, GH = 4L) + (K - 1L)
  if (length(p) < npar) {
    stop('Using ', length(p), ' matching-quantiles to estimate a mixture distribution with ', npar, ' independent parameters. ', 
         'Try increasing `p` (see ?QLMDp for detail).')
  }
  
  qvv <- quantile_vcov(p = p, d = d_obs) 
  qvv_inv <- chol2inv(chol.default(qvv))
  
  parRun <- fmx2dbl(init)
  id_constr <- fmx_constraint_user(distname = distname, K = K, user = constraint)
  # the constraint is from user-specification, instead of from Hoaglin's starting value (i.e. h=0 for negative slope)
  has_constr <- (length(id_constr) > 0L)
  
  # ?stats::optim does not allow Inf in `par`
  par_init <- if (has_constr) parRun[-id_constr] else parRun
  if (any(is.infinite(par_init))) {
    par_init[(par_init < 0) & is.infinite(par_init)] <- -5 
    # .Machine$double.min.exp too small (exp(h) trapped at 0) # exp(-5) close to 0 enough
    # e.g., if h=0 was given in Hoaglin's starting value (due to negative slope), but not as a user-specified constraint,
    # I will use exp(-5)=0.0067 as the starting value, which is passed into ?stats::optim
    if (any((par_init > 0) & is.infinite(par_init))) stop('+Inf parameter indicates ???')
  }
  
  # or use this: https://stackoverflow.com/questions/52552143/how-to-save-the-coefficients-for-each-optim-iteration
  max_return <- .Machine$double.xmax
  
  fn <- if (K == 1L) {
    
    switch(distname, norm = function(x) {
      q <- qnorm(p, mean = x[1L], sd = exp(x[2L]), lower.tail = TRUE, log.p = FALSE)
      if (any(is.infinite(q))) stop('1-comp normal, infinite `q` should not be returned from any `p` between 0 to 1, see `qnorm(.Machine$double.eps)`') # return(max_return)
      return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      
    }, GH = function(x) {
      if (has_constr) parRun[-id_constr] <- x else parRun <- x
      q <- qGH(p, A = parRun[1L], B = exp(parRun[2L]), g = parRun[3L], h = exp(parRun[4L]), lower.tail = TRUE, log.p = FALSE)
      if (any(is.infinite(q))) return(max_return) # stop('interval problem again?') #
      return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      
    })
    
  } else {
    
    Kseq <- seq_len(K)
    Kseq1 <- seq_len(K - 1L)

    switch(distname, norm = {
      id_w <- 2L*K + Kseq1
      function(x) {
        .pM <- array(x[seq_len(2L*K)], dim = c(K, 2L)) # only first 2K elements
        t_w <- t.default(pmlogis_first(x[id_w]))
        sdinv <- 1 / exp(.pM[,2L])
        eff <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L]))) * sdinv
        q <- vuniroot2(y = p, f = function(q) { # essentially \link{pfmx}
          z <- tcrossprod(sdinv, q) - eff
          c(t_w %*% pnorm(z))
        }, interval = interval)
        if (any(is.infinite(q))) return(max_return) # stop('interval problem again?') # 
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }
      
    }, GH = {
      id_w <- 4L*K + Kseq1
      function(x) {
        if (has_constr) parRun[-id_constr] <- x else parRun <- x
        .pM <- array(parRun[seq_len(4L*K)], dim = c(K, 4L)) # only first 4K elements
        t_w <- t.default(pmlogis_first(parRun[id_w]))
        g <- .pM[,3L]
        h <- exp(.pM[,4L])
        sdinv <- 1 / exp(.pM[,2L])
        eff <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L]))) * sdinv
        q <- vuniroot2(y = p, f = function(q) { # essentially \link{pfmx}
          z <- q0 <- tcrossprod(sdinv, q) - eff
          for (i in Kseq) z[i,] <- .qGH2z(q0 = q0[i,], g = g[i], h = h[i], tol = tol, maxiter = maxiter)
          c(t_w %*% pnorm(z))
        }, interval = interval, tol = tol, maxiter = maxiter)
        if (any(is.infinite(q))) return(max_return)
        #if (any(is.infinite(q))) {
        #  print(.pM)
        #  print(p)
        #  print(q)
        #  print(interval)
        #  stop('interval problem again?')
        #}
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }
      
    })
  }
  
  #repeat {
    y <- optim(par = par_init, fn = fn, ...)
    if (FALSE) {
      fn(par_init)
      fn(y$par)
    }
    if (isTRUE(all.equal.numeric(y$par, par_init))) {
      print(init)
      # ret <<- init; x <<- x;
      stop('stats::optim not working (most likely due to poor starting value)')
      # autoplot.fmx(ret, data = x)
      # QLMDe(x, distname = 'GH', K = 2L, init = ret)
    }
    #cat('run stats::optim again..\n')
    # sometimes ?stats::optim will not move.  Weird
  #}
  
  
  if (has_constr) {
    parRun[-id_constr] <- y$par
  } else parRun <- y$par
  tmp <- dbl2fmx(x = parRun, K = K, distname = distname)

  # variance-covariance of internal estimates
  #n <- length(x)
  # `q`: theoretical quantiles
  q <- qfmx(p, distname = distname, K = K, parM = tmp$parM, w = tmp$w, interval = interval, lower.tail = TRUE, log.p = FALSE)
  # `d`: densities at theoretical (\hat{\theta}) quantiles
  d <- dfmx(q, distname = distname, K = K, parM = tmp$parM, w = tmp$w, log = FALSE)
  # {density at \hat{\theta}} and {kernel density} may have different {=0} status.
  tol <- sqrt(sqrt(.Machine$double.eps))
  if (all(d0 <- (abs(d) < tol))) {
    #cat('malformed estimates with all-0 densities\n')
    #return(invisible())
    stop('do not allow this to happen')
  }
  if (any(d0)) {
    p1 <- p[!d0]
    q1 <- q[!d0]
    d1 <- d[!d0]
  } else {
    p1 <- p
    q1 <- q
    d1 <- d
  }
  .meat <- quantile_vcov(p = p1, d = d1) # V_(\hat{theta})
  # in theary, we should use V_{true theta}, but no one knows true theta in practice
  # so I am using V_{\hat_theta}
  # now I want to use V_{empirical} (`@quantile_vv` is a slot of \linkS4class{fmx_QLMDe}), can I?
  q_gr <- qfmx_gr(p = p1, distname = distname, K = K, parM = tmp$parM, w = tmp$w)
  if (!length(q_gr)) {
    int_vv <- array(NA_real_, dim = c(0L, 0L)) # exception handling
  } else {
    .bread <- tryCatch(crossprod_inv(q_gr) %*% t.default(q_gr), error = as.null.default)
    if (!length(.bread)) {
      int_vv <- array(NA_real_, dim = c(0L, 0L)) # exception handling
    } else {
      int_vv <- .bread %*% tcrossprod(.meat, .bread) / length(x) # interval variance-covariance
      if (anyNA(int_vv)) stop('do not allow NA in `int_vv`')
      if (any(diag(int_vv) < 0)) stop('diagonal terms of VarCov < 0 ?!')
      dimnames(int_vv) <- list(dimnames(q_gr)[[2L]], dimnames(q_gr)[[2L]])
    }
  }
  # end of vv
  
  
  ret <- new(
    Class = 'fmx_QLMDe',
    data = x, data.name = data.name,
    distname = distname, parM = tmp$parM, w = tmp$w,
    quantile_vv = qvv,
    vcov = int_vv,
    epdf = x_epdf,
    p = p,
    init = init,
    optim = y
  )
  
  if (!setequal(attr(fmx_constraint(ret), which = 'user', exact = TRUE), constraint)) {
    #message('not handling constrants correctly')
    # indicates ?stats::optim did not move
    #return(invisible())
    stop('should not happen')
  }
  
  return(ret)
  
}



QLMDe_interval <- function(x, extend_interval = 10, ...) {
  xmin <- min(x)
  xmax <- max(x)
  xdiff <- xmax - xmin
  if (length(extend_interval) != 1) stop('...')
  # interval <- c(min(x), max(x)) # old
  c(xmin - extend_interval * xdiff, xmax + extend_interval * xdiff)
}


# 'number of parameters'; this is \strong{not} 'degree-of-freedom'
npar_fmx <- function(x) {
  dm <- dim(x@parM)
  (dm[2L] + 1L) * dm[1L] - 1L - length(attr(fmx_constraint(x), which = 'user', exact = TRUE))
}






# ?base::print
#' @export
print.fmx_QLMDe <- function(x, ...) {
  parM <- x@parM
  K <- dim(parM)[1L]
  dimnames(parM)[[1L]] <- paste0(seq_len(K), '-comp.')
  parM[] <- sprintf(fmt = '%.2f', parM)
  obj <- if (K == 1L) parM else cbind(parM, w = sprintf(fmt = '%.1f%%', x@w*1e2))
  heading <- paste0(K, '-Component Mixture of ', switch(x@distname, norm = 'Normal', GH = 'Tukey\'s G-&-H'), ' Distribution')
  
  ci <- confint.fmx_QLMDe(x, internal = FALSE)
  id_constr <- fmx_constraint(x)
  if (length(ci) && !anyNA(ci)) {
    ci0 <- sprintf(fmt = '(%.2f~%.2f)', ci[,1L], ci[,2L])
    if (length(id_constr)) {
      obj[id_constr] <- '.'
      obj[-id_constr] <- paste(obj[-id_constr], ci0)
    } else obj[] <- paste(obj, ci0)
    heading <- paste0(heading,  ' (w. 95% Confidence Intervals)')
  } else heading <- paste0('Malformed ', heading)
  
  cat('\n ', heading, '\n\n', sep = '')
  print.default(obj, quote = FALSE)
  if (length(id_constr)) cat('\nwhere ', sQuote('.'), ' denotes an enforced constraint\n', sep = '')
  cat('\n')
  
  #cat('\nstats::optim iter.', object@optim$counts['function'], '\n')
  # delta_k = A_k - A_{k-1}
  # no longer print this: if (K > 1L) cat('where \u0394\u2096 = A\u2096 - A\u2096\u208b\u2081, k =', paste0(2:K, collapse = ','), '\n\n') 
  if (inherits(aod <- attr(x, which = 'anova', exact = TRUE), what = 'anova')) {
    print(aod) # ?stats:::print.anova
    cat('\n')
  }
  
  print(autoplot.fmx(x))
  cat('\n')
  return(invisible(x))
}




  

# reparameterization of \linkS4class{fmx} object
# A1 -> A1
# A2 -> A1 + exp(log(d1))
# A_k -> A1 + exp(log(d1)) + .. + exp(log(d_k))
# mixing proportions -> logits
# for 'norm': sd -> log(sd)
# for 'GH': B -> log(B), h -> log(h)
fmx2dbl <- function(x, distname = x@distname, parM = x@parM, K = dim(parM)[1L], w = x@w, ...) { 
  # no longer used in compute intensive algorithms
  w_val <- qmlogis_first(w) # \code{K == 1L} will return \code{numeric(0)}
  if (!all(is.finite(w_val))) stop('NA or Inf in proportion indicated degenerated mixture (one or more component has 0% mixture proportion)')
  w_nm <- if (K == 1L) character() else paste0('logit', 2:K)
  parM[, id] <- log(parM[, (id <- transLog(distname))])
  argnm <- switch(distname, norm = c('mean', 'sdlog'), GH = c('A', 'Blog', 'g', 'hlog'), stop('write more'))
  if (K > 1L) {
    # when some log(d) is negative and has too great absolute value, 
    # exp(log(d)) is numerically 0, and `parM` will not be strictly increasing
    # in the method below we wont be able to retrieve the real log(d), instead get -Inf.
    # we will mark this as constraint in ?qfmx_gr, ?confint.fmx, and ?print.fmx
    parM[2:K, 1L] <- log(parM[2:K, 1L] - parM[1:(K-1L), 1L])
    locnm <- c(paste0(argnm[1L], 1L), paste0('log\u0394', seq_len(K)[-1L]))
  } else locnm <- paste0(argnm[1L], seq_len(K))
  out <- c(parM, w_val)
  names(out) <- c(locnm, paste0(rep(argnm[-1L], each = K), seq_len(K)), w_nm)
  return(out)
}

# invert of ?fmx2dbl
dbl2fmx <- function(x, K, distname, argnm = dist_anm(distname), ...) {
  # not used in ?stats::optim inside \link{QLMDe} (I wrote separate functions for 'GH' and 'norm')
  # only used in \link{qfmx_gr}, which is only used in \link{vcov.fmx_QLMDe}, therefore not compute intensive
  # can use `argnm = NULL` to save time
  nx <- length(x)
  n_dist <- nx - (K - 1L) # K == 1L or not
  w <- if (K == 1L) 1 else unname(pmlogis_first(x[(n_dist + 1L):nx]))
  .pM <- array(x[seq_len(n_dist)], dim = c(K, n_dist/K), dimnames = list(NULL, argnm)) # not compute intensive..
  .pM[,id] <- exp(.pM[, (id <- transLog(distname)), drop = FALSE])
  if (K > 1L) .pM[,1L] <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L])))
  list(parM = .pM, w = w) # much faster than ?base::cbind
}

transLog <- function(distname) {
  switch(distname, norm = {
    2L # `sdlog`
  }, GH = {
    c(2L, 4L) # `Blog` & `hlog`
  }, stop('distribution', sQuote(distname), 'not supported yet'))
}




# ?stats::residuals 
# @export
# residuals.fmx_QLMDe <- function(object, ...) stop('useful?')






#' @title Variance-Covariance of Quantiles
#' 
#' @description
#' 
#' Computes the variance-covariance matrix of quantiles based on Theorem 1 and 2 of Mosteller (1946).
#' 
#' @param p \link[base]{numeric} vector of cumulative probabilities at the given quantiles
#' 
#' @param d \link[base]{numeric} vector of probability densities at the given quantiles
#' 
#' @details 
#' 
#' The end user should make sure no densities too close to 0 is included in argument \code{d}.
#' 
#' \link{quantile_vcov} must not be used in a compute-intensive way.
#' 
#' @return 
#' 
#' The variance-covariance \link[base]{matrix} of quantiles based on Mosteller (1946).
#' 
#' @references 
#' Frederick Mosteller. On Some Useful "Inefficient" Statistics. 
#' \emph{The Annals of Mathematical Statistics}, 17 (4) 377-408, December, 1946. 
#' \doi{10.1214/aoms/1177730881}
#' 
#' @export
quantile_vcov <- function(p, d) {
  # do the check on d=0 in \code{\link{QLMDe}}
  if (anyNA(p) || anyNA(d)) stop('no NA allowed in probability nor density')
  if ((n <- length(p)) != length(d)) stop('p and d must match in length')

  fs <- tcrossprod(d, d) # 'matrix'
  p_c <- array(p, dim = c(n,n)) # 'p on cols'
  p_r <- t.default(p_c) # 'p on rows'
  p_min <- pmin.int(p_r, p_c) # vector!
  p_max <- pmax.int(p_r, p_c)
  vv <- p_min * (1 - p_max) / fs # back to 'matrix'
  return(vv)
}


# gradient of ?qfmx with respect to the \strong{unconstraint} parameters.
# only used in ?vcov.fmx_QLMDe
qfmx_gr <- function(
  dist, # can be missing
  p = stop('must provide `p`'), 
  # constant given `dist`
  distname = dist@distname, parM = dist@parM, K = dim(parM)[1L], w = dist@w,
  interval = qfmx_interval(distname = distname, parM = parM, K = K, w = w, p = c(1e-5, 1-1e-5)),
  ...
) {
  x_skeleton <- fmx2dbl(distname = distname, parM = parM, K = K, w = w)
  has_constr <- (length(id_constr <- fmx_constraint(distname = distname, K = K, parM = parM)) > 0L)
  # ?fmx2dbl may have log(d) being -Inf; see ?fmx2dbl
  x_dbl <- if (has_constr) x_skeleton[-id_constr] else x_skeleton
  x_nm <- names(x_dbl) # stopifnot(all(make.names(x_nm) == x_nm)); # unicode symbols considered legal :)
  if (any(is.infinite(x_dbl))) {
    # stop('wont happen since implementation of constraint') # actually could still happen
    return(invisible())
  }

  if (!is.numeric(interval) || length(interval) != 2L || !all(is.finite(interval))) stop('`interval` must not contain NA nor Inf')

  x_sbl <- lapply(x_nm, FUN = as.symbol)
  qfmx_fn <- as.function.default(c(setNames(rep(x = alist(. = ), times = length(x_nm)), nm = x_nm), as.call(c(
    quote(`{`), 
    quote(x <- x_skeleton),
    if (has_constr) {
      as.call(c(quote(`<-`), quote(x[-id_constr]), as.call(c(quote(c), x_sbl))))
    } else as.call(c(quote(`<-`), quote(x), as.call(c(quote(c), x_sbl)))),
    quote(fx <- dbl2fmx(x = x, distname = distname, K = K, argnm = NULL)),
    quote(qfmx(p = p, parM = fx$parM, w = fx$w, interval = interval, distname = distname, K = K))
  ))))
  
  ret <- tryCatch(expr = {
    # constants `x_skeleton`, `p`, `interval`, `distname`, `K` does not need to be assigned to `rho` :))
    attr(numericDeriv(
      expr = as.call(c(quote(qfmx_fn), x_sbl)), 
      theta = x_nm, 
      central = TRUE, # needed, otherwise could (very rarely) have gradient for `logitw2` all 0
      rho = list2env(as.list.default(x_dbl), parent = environment())), which = 'gradient', exact = TRUE)
  }, error = function(e) { # very rare
    cat('stats::numericDeriv in `qfmx_gr` error.\n')
    #array(NA_real_, dim = c(length(p), length(x_nm)))
    return(invisible())
  })
  if (length(ret)) dimnames(ret)[[2L]] <- x_nm
  return(ret)
}










#' @title Inference for Quantile Least Mahalanobis Distance estimates
#' 
#' @description 
#' 
#' Additional methods of class \linkS4class{fmx} and/or \linkS4class{fmx_QLMDe},
#' for generic functions defined in \pkg{stats} package.
#' 
#' @param object an \linkS4class{fmx} or \linkS4class{fmx_QLMDe} object
#' 
#' @param data \link[base]{double} vector, actual observations
#' 
#' @param internal \link[base]{logical} scalar, either for the user-friendly parameters (\code{FALSE}, default)
#' (e.g., \code{mean,sd} for normal mixture, and \code{A,B,g,h} for Tukey's \eqn{g}-and-\eqn{h} mixture), or
#' for the internal/unconstrained parameters (\code{TRUE}).
#' 
#' @param level confidence level, default \eqn{95\%}.
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @details 
#' 
#' The inference for the user-friendly parameters is obtained via delta-method. 
#' 
#' @return 
#' 
#' \link{nobs.fmx_QLMDe} returns an \link[base]{integer} scalar indicating the sample size of
#' the observations used in \link{QLMDe} estimation. 
#' 
#' \link{logLik.fmx} returns a \link[stats]{logLik} object indicating the log-likelihood.
#' An additional attribute \code{attr(, 'logl')} indicates the point-wise log-likelihood, 
#' to be use in Vuong's closeness test.
#' 
#' \link{coef.fmx} returns the estimates of the user-friendly parameters (\code{parm = 'user'}), 
#' or the internal/unconstrained parameters (\code{parm = 'internal'}).
#' 
#' \link{vcov.fmx_QLMDe} returns 
#' the approximate asymptotic variance-covariance matrix of the user-friendly parameters via delta-method (\code{parm = 'user'}), 
#' or the asymptotic variance-covariance matrix of the internal/unconstrained parameters (\code{parm = 'internal'}). 
#' 
#' \link{confint.fmx_QLMDe} returns the Wald-type confidence intervals based on the user-friendly parameters (\code{parm = 'user'}),
#'  or the internal/unconstrained parameters (\code{parm = 'internal'}).
#' 
#' When the distribution has constraints on one or more parameters, none of \code{\link{coef.fmx}}, \link{vcov.fmx_QLMDe} and 
#' \link{confint.fmx_QLMDe} return the corresponding values only for the constrained parameters.
#' 
#' @seealso \link[stats]{nobs} \link[stats]{logLik} \link[stats]{coef} \link[stats]{vcov}
#' \link[stats]{confint}
#' 
#' @name S3_fmx_QLMDe
#' @export
vcov.fmx_QLMDe <- function(object, internal = FALSE, ...) {
  
  int_vv <- object@vcov
  if (internal) return(int_vv)
  
  if (!length(int_vv)) return(int_vv) # wont be able to computer user-vcov if internal-vcov is wrong
  
  distname <- object@distname
  parM <- object@parM
  K <- dim(parM)[1L]
  int_nm <- dimnames(int_vv)[[1L]]
  
  int_p <- fmx2dbl(object) # internal parameters
  anm <- dist_anm(distname)
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



#' @rdname S3_fmx_QLMDe
#' @export
coef.fmx <- function(object, internal = FALSE, ...) {
  anm <- dist_anm(object@distname)
  K <- dim(object@parM)[1L]
  cf0 <- if (internal) fmx2dbl(object) else if (K == 1L) {
    setNames(c(object@parM), nm = c(t.default(outer(anm, 1:K, FUN = paste0))))
  } else {
    setNames(c(object@parM, object@w), nm = c(t.default(outer(c(anm, 'w'), 1:K, FUN = paste0))))
  }
  if (!length(id_constr <- fmx_constraint(object))) return(cf0)
  return(cf0[-id_constr])
}


#' @rdname S3_fmx_QLMDe
#' @export
confint.fmx_QLMDe <- function(object, ..., level = .95) {
  # essentially ?stats::confint.default
  cf <- coef.fmx(object, ...)
  if (!length(vv <- vcov.fmx_QLMDe(object, ...))) return(invisible())
  ses <- sqrt(diag(vv))
  p1 <- (1 - level) / 2
  p <- c(p1, 1 - p1)
  ret <- cf + ses %*% t.default(qnorm(p))
  dimnames(ret) <- list(names(cf), sprintf('%.1f%%', 1e2*p))
  return(ret)
}



#' @rdname S3_fmx_QLMDe
#' @export
logLik.fmx <- function(object, data = object@data, ...) {
  # for developer to batch-calculate AIC/BIC quickly
  if (inherits(object, 'fmx_QLMDe') && (nobjF <- length(objF <- attr(object, which = 'objF', exact = TRUE)))) {
    if (inherits(objF[[nobjF]], what = 'logLik')) return(objF[[nobjF]])
  }
  
  logd <- dfmx(x = data, dist = object, log = TRUE, ...)
  if (!all(is.finite(logd))) {
    #print(logd)
    #object <<- object
    #stop('malformed fit (?.dGH has been well debug-ged)')
    # very likely to be `B = 0`
    # do not stop.  settle with -Inf log-likelihood
  }
  out <- sum(logd)
  attr(out, which = 'logl') <- logd # additional attributes; needed in Vuong's test
  attr(out, which = 'nobs') <- length(data)
  
  # https://en.wikipedia.org/wiki/Akaike_information_criterion
  # ?stats::logLik
  # ?stats:::AIC.default
  # attr(, 'df') is the number of (estimated) parameters in the model.
  attr(out, which = 'df') <- npar_fmx(object)
  
  class(out) <- 'logLik'
  return(out)
}

# ?stats:::AIC.logLik





#' @rdname S3_fmx_QLMDe
#' @export
nobs.fmx_QLMDe <- function(object, ...) length(object@data)




