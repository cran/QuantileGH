
#' @title Quantile Least Mahalanobis Distance estimates
#' 
#' @description 
#' 
#' The quantile least Mahalanobis distance algorithm (\code{\link{QLMDe}}) estimates the parameters of 
#' single-component or finite mixture distributions   
#' by minimizing the Mahalanobis distance between the vectors of sample and theoretical quantiles.
#' 
#' {Add a paragraph about quantiles}
#' 
#' The initial values of finite mixture distribution are based on trimmed \eqn{k}-means 
#' clustering with re-assignment provided by function \code{\link{clust_fmx}}.
#' 
#' @param x vector of 'numeric' observations
#' 
#' @param distname 'character', name of mixture distribution to be fitted.  Currently supports \code{'norm'} and \code{'GH'}
#' 
#' @param K 'integer', number of components (must use \code{2L} instead of \code{2})
#' 
#' @param data.name 'character', user-assigned name for the observations. Default is the model call of the input observations.
#' 
#' @param p vector of \code{'numeric'} percentiles, at where the sample and theoretical quantiles are to be matched.
#' Default \code{QLMDp(obs = x)}, see \code{\link{QLMDp}} for details.
#' 
#' @param init \code{'fmx'} object, initial values of the optimization algorithm. Default value is provided by function \code{\link{clust_fmx}}.
#' 
#' @param constraint the parameters (\eqn{g} and/or \eqn{h} for Tukey's \eqn{g}-&-\eqn{h} mixture) to be set at 0.  
#' See \code{\link{fmx_constraint}} for details.
#' 
#' @param method,control see \code{\link[stats]{optim}}.
#' 
#' @param ... only potential parameters of \code{\link[stats]{optim}} are allowed.
#' 
#' @details
#' 
#' Quantile Least Mahalanobis Distance estimator (\code{\link{QLMDe}}) fits a single-component or finite mixture distribution 
#' by minimizing the Mahalanobis distance between
#' the theoretical and observed quantiles,
#' using the empirical quantile variance-covariance matrix \code{\link{quantile_vcov}}.
#' 
#' @return An \code{'fmx_QLMDe'} object (see \code{\link{fmx_QLMDe-class}}).
#' 
#' @seealso \code{\link[stats]{optim}}
#' 
#' @examples 
#' 
#' set.seed(1623); (y0n <- QLMDe(rnorm(1e3L), distname = 'norm', K = 2L))
#' \donttest{
#' # slow, masked for package check
#' StepK_fmx(y0n, by = 'logLik') # one-component
#' }
#' 
#' set.seed(1623); x0gh = rGH(1e3L, g = .2, h = .1)
#' \donttest{
#' print(y0gh <- QLMDe(x0gh, distname = 'GH', K = 2L))
#' print(StepK_fmx(y0gh, by = 'logLik')) # really good!
#' }
#' 
#' (d1 <- fmx('norm', mean = c(0, 1.5), sd = .5, w = c(.4, .6)))
#' set.seed(100); hist(x1 <- rfmx(n = 1e3L, dist = d1))
#' \donttest{
#' system.time(tmp <- QLMDe(x1, distname = 'norm', K = 1L))
#' tmp
#' QLMDe(x1, distname = 'norm', K = 2L)
#' y1n <- StepK_fmx(tmp, by = 'logLik')
#' slotNames(y1n)
#' logLik(y1n)
#' AIC(y1n)
#' BIC(y1n)
#' }
#'   
#' (d2 = fmx('GH', A = c(0,6), B = 1.2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(3123); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' \donttest{
#' system.time(y2gh <- QLMDe(x2, distname = 'GH', K = 2L)) # 8~9 secs
#' y2gh@@optim$counts # ~840 iteration
#' y2gh
#' Step_fmx(y2gh, by = 'logLik') # not always identifying true constraint
#' 
#' system.time(y2gh_0 <- QLMDe(x2, distname = 'GH', K = 2L, constraint = c('g1', 'h2')))
#' y2gh_0
#' system.time(y2n <- QLMDe(x2, distname = 'norm', K = 2L))
#' y2n
#' print(StepK_fmx(y2gh_0))
#' stopifnot(AIC(y2gh_0) < AIC(y2n), BIC(y2gh_0) < BIC(y2n))
#' }
#' 
#' @name QLMDe
#' @export
QLMDe <- function(
  x, distname, K, data.name = deparse1(substitute(x)),
  constraint = character(),
  p = QLMDp(obs = x),
  init = clust_fmx(x, distname = distname, K = K, constraint = constraint),
  method = 'Nelder-Mead', 
  control = list(maxit = 1e4L, reltol = 1e-6),
  ...
) {
  
  if (!is.character(distname) || length(distname) != 1L || anyNA(distname) || !nzchar(distname)) stop('distname should be len-1 character')
  if (!is.integer(K) || length(K) != 1L || anyNA(K) || K <= 0L) stop('number of component must be len-1 positive integer')
  
  if (!is.vector(x, mode = 'double')) stop('x must be double vector')
  if (anyNA(x)) stop('do not allow NA in observations')
  if (!is.character(data.name) || length(data.name) != 1L || anyNA(data.name) || !nzchar(data.name)) stop('illegal data.name')
  
  xmin <- min(x)
  xmax <- max(x)
  q_lim <- c(xmin, xmax) # do not think about extending this interval!!!
  method <- match.arg(method, choices = eval(formals(optim)$method))
  if (method == 'Brent') stop('Brent method is only for 1-dim problem')
  
  #p <- QLMDp(obs = x, ...)
  if (anyNA(p)) stop('should not happen')
  npar <- K * switch(distname, norm = 2L, GH = 4L, stop(sQuote(distname), ' not supported')) + (K - 1L)
  
  q_init <- qfmx(p = p, interval = q_lim, distname = distname, K = K, parM = init@parM, w = init@w)
  if (!all(id1 <- is.finite(q_init))) {
    if (!any(id1)) stop('starting values too far away?')
    p <- p[id1]
  }
  q_obs <- quantile(x, probs = p) # observed quantiles, constant in ?stats::optim
  
  x_kern <- density.default(x) # ?stats::approx inside ?stats::density.default
  x_epdf <- approxfun(x = x_kern$x, y = x_kern$y) # another 'layer' of ?stats::approxfun
  # Tingting is not sure whether ?stats::approx \strong{and} ?stats::approxfun will make `x_epdf(q_obs) = 0` more likely
  d_obs <- x_epdf(q_obs) # observed density evaluated at `q_obs`
  tol <- sqrt(sqrt(.Machine$double.eps))
  if (anyNA(d_obs)) stop('should not happen')
  if (all(d0 <- (abs(d_obs) < tol))) stop('must have at least one positive density') # `d_obs` should always be positive, but ?base::abs should not hurt
  if (any(d0)) {
    p <- p[!d0] # quantiles where `density = 0` are excluded from the calculation 
    d_obs <- d_obs[!d0]
    q_obs <- q_obs[!d0]
  }
  
  if (length(p) < npar) {
    stop('Using ', length(p), ' matching-quantiles to estimate a mixture distribution with ', npar, ' independent parameters. ', 
         'Try increasing `p.N` or use more `p.extra` (see ?QLMDp for detail).')
  }
  
  qvv <- quantile_vcov(p = p, d = d_obs) 
  qvv_inv <- chol2inv(chol(qvv))
  
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
    # I will use exp(-5)=0.0067 as the starting value that's being passed into ?stats::optim
    if (any((par_init > 0) & is.infinite(par_init))) stop('Inf parameter indicates ???')
  }
  
  # or use this: https://stackoverflow.com/questions/52552143/how-to-save-the-coefficients-for-each-optim-iteration
  #parTrace <- .Internal(array(par_init, c(1L, length(par_init)), list(NULL, names(par_init)))) # I don't know how to use 'trace' in `control` of \code{stats::optim}
  #env <- environment()
  #optim_ctrl <- list(maxit = maxit, reltol = reltol, parscale = parscale)
  max_return <- .Machine$double.xmax
  
  y <- if (K == 1L) {
    
    if (distname == 'norm') {
      optim(par = par_init, fn = \(x) {
        q <- qnorm(p, mean = x[1L], sd = exp(x[2L]), lower.tail = TRUE, log.p = FALSE)
        if (any(is.infinite(q))) return(max_return)
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }, method = method, control = control, ...) #, hessian = TRUE
    } else if (distname == 'GH') {
      optim(par = par_init, fn = \(x) {
        if (has_constr) parRun[-id_constr] <- x else parRun <- x
        q <- qGH(p, A = parRun[1L], B = exp(parRun[2L]), g = parRun[3L], h = exp(parRun[4L]), lower.tail = TRUE, log.p = FALSE)
        if (any(is.infinite(q))) return(max_return)
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }, method = method, control = control, ...) #, hessian = TRUE
    } else stop('distribution ', sQuote(distname), 'not ready')
    
  } else {
    
    Kseq <- seq_len(K)
    Kseq1 <- seq_len(K - 1L)
    # tol <- .Machine$double.eps^.25 # default of ?rstpm2::vuniroot
    
    # vuniroot_interval <- t.default(array(q_lim, dim = c(2L, length(p))))
    
    if (distname == 'norm') {
      id_w <- 2L*K + Kseq1
      optim(par = par_init, fn = \(x) {
        .pM <- array(x[seq_len(2L*K)], dim = c(K, 2L)) # only first 2K elements
        t_w <- t.default(pmlogis_first(x[id_w]))
        sdinv <- 1 / exp(.pM[,2L])
        eff <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L]))) * sdinv
        #q <- tryCatch(vuniroot(f = \(q) { # essentially \code{\link{pfmx}}
        #  z <- tcrossprod(sdinv, q) - eff # remove .Internal for CRAN
        #  ret <- c(t_w %*% pnorm(z)) 
        #  ret - p # needed for ?rstpm2::vuniroot
        #}, interval = vuniroot_interval)$root, error = identity)
        q <- tryCatch(vuniroot2(y = p, f = \(q) { # essentially \code{\link{pfmx}}
          z <- tcrossprod(sdinv, q) - eff # remove .Internal for CRAN
          ret <- c(t_w %*% pnorm(z)) 
          #ret - p # needed for ?rstpm2::vuniroot
        }, interval = q_lim), error = identity)
        if (inherits(q, what = 'error')) return(max_return)
        return(mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv))
      }, method = method, control = control, ...) #, hessian = TRUE
      
    } else if (distname == 'GH') {
      id_w <- 4L*K + Kseq1
      optim(par = par_init, fn = function(x) {
        #assign('parTrace', value = rbind(parTrace, x), envir = env)
        if (has_constr) parRun[-id_constr] <- x else parRun <- x
        .pM <- array(parRun[seq_len(4L*K)], dim = c(K, 4L)) # only first 4K elements
        t_w <- t.default(pmlogis_first(parRun[id_w]))
        g <- .pM[,3L]
        h <- exp(.pM[,4L])
        sdinv <- 1 / exp(.pM[,2L])
        eff <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L]))) * sdinv
        #q <- tryCatch(vuniroot(f = \(q) { # essentially \code{pfmx}
        #  z <- q0 <- tcrossprod(sdinv, q) - eff
        #  for (i in Kseq) z[i,] <- .qGH2z(q0 = q0[i,], g = g[i], h = h[i])
        #  ret <- c(t_w %*% pnorm(z))
        #  ret - p
        #}, interval = vuniroot_interval)$root, error = identity)
        q <- tryCatch(vuniroot2(y = p, f = \(q) { # essentially \code{pfmx}
          z <- q0 <- tcrossprod(sdinv, q) - eff
          for (i in Kseq) z[i,] <- .qGH2z(q0 = q0[i,], g = g[i], h = h[i])
          ret <- c(t_w %*% pnorm(z))
          #ret - p
        }, interval = q_lim), error = identity)
        if (inherits(q, what = 'error')) return(max_return)
        out_Maha <- mahalanobis_int(x = q, center = q_obs, invcov = qvv_inv)
        if (is.infinite(out_Maha)) stop('should not happen')
        return(out_Maha)
      }, method = method, control = control, ...) #, hessian = TRUE
      
    } else stop('distribution not ready: ', sQuote(distname))
  }
  
  # information of ?stats::optim (R 2021-02-16)
  y$method <- method
  #y$parTrace <- parTrace
  
  if (has_constr) parRun[-id_constr] <- y$par else parRun <- y$par
  ret <- dbl2fmx(x = parRun, K = K, distname = distname)

  new(
    Class = 'fmx_QLMDe',
    data = x, data.name = data.name,
    distname = distname, parM = ret$parM, w = ret$w,
    quantile_vv = qvv,
    epdf = x_epdf,
    p = p,
    init = init,
    optim = y
  )
  
}


# 'number of parameters'; this is \strong{not} 'degree-of-freedom'
npar_fmx <- function(x, ...) {
  dm <- dim(x@parM)
  (dm[2L] + 1L) * dm[1L] - 1L - length(attr(fmx_constraint(x), which = 'user', exact = TRUE))
}


# extending ?stats::logLik
#' @export
logLik.fmx_QLMDe <- function(object, ...) {
  # for developer to batch-calculate AIC/BIC quickly
  #if (inherits(atr <- attr(object, which = 'logLik', exact = TRUE), what = 'logLik')) return(atr) 
  names(attributes(object))
  if (nobjF <- length(objF <- attr(object, which = 'objF', exact = TRUE))) {
    if (inherits(objF[[nobjF]], what = 'logLik')) return(objF[[nobjF]])
  }
  
  parM <- object@parM
  K <- dim(parM)[1L]
  logd <- dfmx(x = object@data, log = TRUE, distname = object@distname, K = K, parM = parM, w = object@w, ...)
  if (anyNA(logd)) stop('mal-shaped fit (?.dGH has been well debug-ged)')
  out <- sum(logd)
  attr(out, 'logl') <- logd # addl. attr. for 'logLik' objects (needed in Vuong's test)
  attr(out, 'nobs') <- length(object@data)
  attr(out, 'npar') <- npar_fmx(object)
  attr(out, 'df') <- attr(out, which = 'nobs', exact = TRUE) - attr(out, which = 'npar', exact = TRUE)
  class(out) <- 'logLik'
  return(out)
}

# ?stats:::AIC.logLik




# ?stats::nobs
#' @export
nobs.fmx_QLMDe <- function(object, ...) length(object@data)












# ?base::print
#' @export
print.fmx_QLMDe <- function(x, ...) {
  parM <- x@parM
  K <- dim(parM)[1L]
  parM[] <- sprintf(fmt = '%.2f', parM)
  obj <- if (K == 1L) parM else cbind(parM, w = sprintf(fmt = '%.1f%%', x@w*1e2))
  heading <- paste0(K, '-Component Mixture of ', x@distname, ' Distribution')
  
  ci <- confint.fmx_QLMDe(x, parm = 'user')
  if (length(ci)) {
    ci0 <- sprintf(fmt = '%.2f~%.2f', ci[,1L], ci[,2L])
    if (length(id_constr <- fmx_constraint(x))) {
      obj[id_constr] <- paste0(obj[id_constr], ' (constraint)')
      obj[-id_constr] <- paste0(obj[-id_constr], ' (', ci0, ')')
    } else obj[] <- paste0(obj, ' (', ci0, ')')
    heading <- paste0(heading,  ' (w. 95% Confidence Intervals)')
  } else heading <- paste0('Malformed ', heading)
  
  dimnames(obj)[[1L]] <- paste0(seq_len(K), '-comp.')
  cat('\n ', heading, '\n\n', sep = '')
  print.default(obj, quote = FALSE)
  cat('\n')
  
  #cat('\nstats::optim iter.', object@optim$counts['function'], '\n')
  # delta_k = A_k - A_{k-1}
  # no longer print this: if (K > 1L) cat('where \u0394\u2096 = A\u2096 - A\u2096\u208b\u2081, k =', paste0(2:K, collapse = ','), '\n\n') 
  if (inherits(aod <- attr(x, which = 'anova', exact = TRUE), what = 'anova')) {
    print(aod) # ?stats:::print.anova
    cat('\n')
  }
  
  print(autoplot.fmx_QLMDe(x))
  cat('\n')
  return(invisible(x))
}


#' @export
`[.fmx_QLMDe` <- function(x, i) stop('do not allow subsetting \'fmx_QLMDe\' object!')







# reparameterization of 'fmx' object
# A1 -> A1
# A2 -> A1 + exp(log(d1))
# A_k -> A1 + exp(log(d1)) + ... + exp(log(d_k))
# mixing proportions -> logits
# for 'norm': sd -> log(sd)
# for 'GH': B -> log(B), h -> log(h)
fmx2dbl <- function(x, distname = x@distname, parM = x@parM, K = dim(parM)[1L], w = x@w, ...) { 
  # no longer used in compute intensive algorithms
  # logit(proportion)
  # proper log-transformed arguments
  w_val <- qmlogis_first(w) # \code{K == 1L} will return \code{numeric(0)}
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
  # not used in ?stats::optim inside \code{\link{QLMDe}} (I wrote separate functions for 'GH' and 'norm')
  # only used in \code{\link{qfmx_gr}}, which is only used in \code{\link{vcov.fmx_QLMDe}}, therefore not compute intensive
  nx <- length(x)
  if (K == 1L) {
    n_dist <- nx
    w <- 1
  } else {
    n_logits <- K - 1L
    n_dist <- nx - n_logits
    w <- unname(pmlogis_first(x[(n_dist + 1L):nx]))
    if (anyNA(w)) {
      print(x[(n_dist + 1L):nx])
      stop('these multinomial logits creates NA proportions?')
    }
  }
  # .pM <- .Internal(array(x, c(K, n_dist/K), NULL)) # extra elements dropped!!
  .pM <- array(x[seq_len(n_dist)], dim = c(K, n_dist/K)) # not compute intensive..
  .pM[,id] <- exp(.pM[, (id <- transLog(distname))])
  if (K > 1L) .pM[,1L] <- cumsum(c(.pM[1L,1L], exp(.pM[2:K,1L])))
  dimnames(.pM)[[2L]] <- argnm
  list(parM = .pM, w = w) # much faster than ?base::cbind
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
#' @param p,d probability and density at the given quantiles
#' 
#' @details 
#' 
#' The end user should make sure no densities too close to 0 is included in argument \code{d}.
#' 
#' @return 
#' 
#' The variance-covariance \code{'matrix'} of quantiles based on Mosteller (1946).
#' 
#' @references 
#' Frederick Mosteller. On Some Useful "Inefficient" Statistics. The Annals of Mathematical Statistics, 17 (4) 377-408, December, 1946. 
#' \doi{10.1214/aoms/1177730881}
#' 
#' @export
quantile_vcov <- function(p, d) {
  
  # stop('do the check on d=0 in QLMDe now')
  if (anyNA(p) || anyNA(d)) stop('no NA allowed in probability nor density')
  # if (any(duplicated_allequal(p))) stop('numerically duplicated p; programming error')
  if ((n <- length(p)) != length(d)) stop('p and d must match in length')
  #if (all(d0 <- (abs(d) < tol))) stop('must have at least one positive density')
  #p <- p[!d0] # quantiles where `density = 0` are excluded from the calculation 
  #d <- d[!d0]
  
  # not compute-intensive
  fs <- tcrossprod(d, d) # 'matrix'
  p_c <- array(p, dim = c(n,n)) # 'p on cols'
  p_r <- t.default(p_c) # 'p on rows'
  # end of not compute-intensive
  p_min <- pmin.int(p_r, p_c) # vector!
  p_max <- pmax.int(p_r, p_c)
  vv <- p_min * (1 - p_max) / fs # back to 'matrix'
  # attr(vv, which = 'd0') <- d0
  return(vv)
}


# gradient of ?qfmx with respect to \strong{unconstraint} parameters.
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
  if (any(is.infinite(x_dbl))) stop('wont happen since implementation of constraint')

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
  
  # constants `x_skeleton`, `p`, `interval`, `distname`, `K` does not need to be assigned to `rho` :))
  out <- numericDeriv(
    expr = as.call(c(quote(qfmx_fn), x_sbl)), 
    theta = x_nm, 
    rho = list2env(as.list.default(x_dbl), parent = environment())
  ) # tryCatch(, error = identity); let err
  #if (inherits(out, 'error')) return(invisible()) # very rare
  out_gr <- attr(out, which = 'gradient', exact = TRUE)
  dimnames(out_gr)[[2L]] <- x_nm
  return(out_gr)
}










#' @title Inference for Quantile Least Mahalanobis Distance estimates
#' 
#' @description S3 method for \code{\link[stats]{vcov}}, \code{\link[stats]{coef}} and \code{\link[stats]{confint}}.
#' 
#' @param object an \code{'fmx'} or \code{'fmx_QLMDe'} object
#' 
#' @param parm \code{'character'} value. Use \code{'user'} for the user-friendly parameters (e.g., \code{mean,sd} for normal, and \code{A,B,g,h} for Tukey's \eqn{g}-and-\eqn{h}), via delta-method.
#' Use \code{'internal'} for the internal/unconstrained parameters.
#' 
#' @param level confidence level, default \eqn{95\%}.
#' 
#' @param ... place holder for S3 naming convention
#' 
#' @return 
#' 
#' \code{\link{coef.fmx}} returns the estimates of the user-friendly parameters (\code{parm = 'user'}), 
#' or the internal/unconstrained parameters (\code{parm = 'internal'}).
#' 
#' \code{\link{vcov.fmx_QLMDe}} returns 
#' the approximate asymptotic variance-covariance matrix of the user-friendly parameters via delta-method (\code{parm = 'user'}), 
#' or the asymptotic variance-covariance matrix of the internal/unconstrained parameters (\code{parm = 'internal'}). 
#' 
#' \code{\link{confint.fmx_QLMDe}} returns the Wald-type confidence intervals based on the user-friendly parameters (\code{parm = 'user'}),
#'  or the internal/unconstrained parameters (\code{parm = 'internal'}).
#' 
#' When the distribution has constraints on one or more parameters, none of \code{\link{coef.fmx}}, \code{\link{vcov.fmx_QLMDe}} and 
#' \code{\link{confint.fmx_QLMDe}} will return the corresponding values only for the constrained parameters.
#' 
#' @name S3_fmx_QLMDe
#' @method vcov fmx_QLMDe
#' @export
vcov.fmx_QLMDe <- function(object, parm = c('user', 'internal'), ...) {
  # not compute-intensive
  parm <- match.arg(parm)
  distname <- object@distname
  parM <- object@parM
  K <- dim(parM)[1L]
  w <- object@w
  x <- object@data
  n <- length(x)
  interval <- c(min(x), max(x))
  p <- object@p
  
  # `q`: theoretical quantiles
  q <- qfmx(p, distname = distname, K = K, parM = parM, w = w, interval = interval, lower.tail = TRUE, log.p = FALSE)
  # if (any(xid <- is.infinite(q))) then allow infinite
  
  # `d`: densities at theoretical quantiles
  d <- dfmx(q, distname = distname, K = K, parM = parM, w = w, log = FALSE) # density at \hat{\theta}
  # {density at \hat{\theta}} and {kernel density} may have different {=0} status.
  # stop('now quantile_vv is a slot of fmx_QLMDe')
  .meat <- quantile_vcov(p = p, d = d) # V_(\hat{theta})
  # in theary, we should use V_{true theta}, but no one knows the true theta in practice
  # so I am using V_{\hat_theta}
  # now I want to use V_{empirical}, can I?
  
  q_gr <- tryCatch(qfmx_gr(p = p, distname = distname, K = K, parM = parM, w = w), error = as.null.default)
  # in ?qfmx_gr, `interval` could have error, although very rarely. 
  if (!length(q_gr)) return(array(NA_real_, dim = c(0L, 0L))) 
  .bread <- crossprod_inv(q_gr) %*% t.default(q_gr)
  int_vv <- .bread %*% tcrossprod(.meat, .bread) / n
  if (any(diag(int_vv) < 0)) stop('diagonal terms of VarCov < 0 ?!')
  dimnames(int_vv) <- list(int_nm <- dimnames(q_gr)[[2L]], int_nm)
  if (parm == 'internal') return(int_vv)
  
  int_p <- fmx2dbl(object)
  anm <- dist_anm(distname)
  n_anm <- length(anm)
  user_nm <- c(t.default(outer(c(anm, if (K > 1L) 'w'), 1:K, FUN = paste0)))
  jacob <- array(0, dim = c(length(int_p), length(user_nm)), dimnames = list(names(int_p), user_nm))
  # location parameters A_1, A_2, ..., A_k -> A_1, d_2, ..., d_k
  jacob[1L, 1:K] <- 1
  if (K > 1L) for (k in 2:K) jacob[k, k:K] <- exp(int_p[k])
  # end of location parameters
  if (K > 1L) {
    # mixture parameters w_1, ..., w_k -> pi_2, ..., pi_k
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
coef.fmx <- function(object, parm = c('user', 'internal'), ...) {
  anm <- dist_anm(object@distname)
  K <- dim(object@parM)[1L]
  cf0 <- switch(match.arg(parm), internal = fmx2dbl(object), user = {
    if (K > 1L) {
      setNames(c(object@parM, object@w), nm = c(t.default(outer(c(anm, 'w'), 1:K, FUN = paste0))))
    } else setNames(c(object@parM), nm = c(t.default(outer(anm, 1:K, FUN = paste0))))
  })
  if (!length(id_constr <- fmx_constraint(object))) return(cf0)
  return(cf0[-id_constr])
}


#' @rdname S3_fmx_QLMDe
#' @export
confint.fmx_QLMDe <- function(object, parm = c('internal', 'user'), level = .95, ...) {
  parm <- match.arg(parm)
  return(confint_int(cf = coef.fmx(object, parm = parm), vv = vcov.fmx_QLMDe(object, parm = parm), level = level))
}


# @rdname show_S4
# @export
# setMethod(show, signature(object = 'fmx_QLMDe'), definition = \(object) print.fmx(object))
# inherits from `show.fmx`







