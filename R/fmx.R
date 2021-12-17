

#' @title The Finite Mixture Distribution
#' 
#' @description 
#' 
#' Density function, distribution function, quantile function and random generation for a finite mixture distribution 
#' with normal or Tukey's \eqn{g}-&-\eqn{h} components.
#' 
#' @param x,q vector of quantiles.
#' 
#' @param p vector of probabilities.
#' 
#' @param n number of observations.
#' 
#' @param dist \code{'fmx'} object, representing a finite mixture distribution
#' 
#' @param log,log.p logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{log(p)}.
#' 
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{Pr(X\le x)} otherwise, \eqn{Pr(X>x)}.
#' 
#' @param interval interval for root finding (see \code{\link[rstpm2]{vuniroot}})
#' 
#' @param distname,K,parM,w auxiliary parameters, whose default values are determined by argument \code{dist}, 
#' see \code{\link{fmx-class}} for details.
#' The user-specified vector of \code{w} does not need to sum up to 1; \code{w/sum(w)} will be used internally.
#' 
#' @param ... in \code{\link{fmx}}, these are the vectors of distribution parameters
#' 
#' @details 
#' 
#' \code{\link{fmx}} creates an S4 \code{'fmx'} object which specifies the parameters of a finite mixture distribution.
#' See \code{\link{fmx-class}} for details.
#' 
#' \code{\link{dfmx}} gives the density, which is required in plotting and calculating
#' the log-likelihood. 
#' 
#' \code{\link{pfmx}} gives the distribution function.
#' 
#' \code{\link{qfmx}} gives the quantile function, by solving \code{\link{pfmx}} by \code{\link[rstpm2]{vuniroot}}.
#' One major challenge when dealing with the finite mixture of Tukey's \eqn{g}-&-\eqn{h} family distribution
#' is that Brent–Dekker's method needs to be performed in both \code{\link{pGH}} and \code{\link{qfmx}}, 
#' i.e. 'two layers' of root-finding algorithm.
#' 
#' \code{\link{rfmx}} generates random deviates. 
#' 
#' @return 
#' 
#' \code{\link{fmx}} returns an S4 \code{'fmx'} object which specifies the parameters of a finite mixture distribution.
#' See \code{\link{fmx-class}} for details.
#' 
#' \code{\link{dfmx}} returns a vector of probability density values of an S4 \code{'fmx'} object at specified quantiles \code{x}.
#' 
#' \code{\link{pfmx}} returns a vector of cumulative probability values of an S4 \code{'fmx'} object at specified quantiles \code{q}.
#' 
#' \code{\link{qfmx}} returns a vector of quantiles of an S4 \code{'fmx'} object, based on specified cumulative probabilities \code{p}.
#' 
#' \code{\link{rfmx}} generates random deviates of an S4 \code{'fmx'} object.
#' 
#' @examples 
#' 
#' x = (-3):7
#' 
#' (e1 = fmx('norm', mean = c(0,3), sd = c(1,1.3), w = c(1, 1)))
#' # 'fmx' object `e1`, a two-component 50%-50% mixture of normal
#' isS4(e1) # 'fmx' is S4 object
#' slotNames(e1) # 'slot' is similar to S3 object component 
#' e1@@parM # component parameters
#' plot(e1) # using vanilla R
#' autoplot(e1) # using ggplot2 package
#' hist(rfmx(n = 1e3L, dist = e1), main = '1000 obs from e1')
#' # generate a sample of size 1e3L from mixture distribution `e1`
#' round(dfmx(x, dist = e1), digits = 3L)
#' round(p1 <- pfmx(x, dist = e1), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p1, dist = e1), x, tol = 1e-4))
#' 
#' (e2 = fmx('GH', A = c(0,3), g = c(.2, .3), h = c(.2, .1), w = c(2, 3)))
#' # 'fmx' object `e2`, a two-component 40%-60% mixture of Tukey
#' hist(rfmx(n = 1e3L, dist = e2), main = '1000 obs from e2') 
#' round(dfmx(x, dist = e2), digits = 3L)
#' round(p2 <- pfmx(x, dist = e2), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p2, dist = e2), x, tol = 1e-4))
#' 
#' (e3 = fmx('GH', A = 0, g = .2, h = .2)) # one-component Tukey
#' hist(rfmx(1e3L, dist = e3))
#' hist(rGH(n = 1e3L, A = 0, g = .2, h = .2))
#' # identical (up to random seed); but ?rfmx has much cleaner code
#' round(dfmx(x, dist = e3), digits = 3L)
#' round(p3 <- pfmx(x, dist = e3), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p3, dist = e3), x, tol = 1e-4))
#' 
#' @name fmx
#' @export
fmx <- function(distname, w = 1, ...) {
  if (!is.character(distname) || (length(distname) != 1L) || anyNA(distname) || !nzchar(distname)) stop('distname must be len-1 char')
  anm <- dist_anm(distname)
  arg <- list(...)[anm]
  if (!any(la <- lengths(arg, use.names = FALSE))) {
    cat('no user-provided args, use default arguments of', sQuote(paste0('d', distname)), '\n')
  }
  if (length(la_m <- la[la > 1L])) { # multiple component
    if (!all(duplicated.default(la_m)[-1L])) stop('user-provided args must be of same length')
    if (any(id_la1 <- (la == 1L))) arg[id_la1] <- lapply(arg[id_la1], FUN = rep, times = la_m[1L])
  }
  if (any(id_fa <- !la)) {
    farg <- formals(paste0('d', distname))[anm[id_fa]]
    if (any(id_fa0 <- !vapply(farg, FUN = nzchar, FUN.VALUE = NA))) stop('formal density argument not available: ', sQuote(anm[id_fa0]))
    arg[id_fa] <- if (!length(la_m)) farg else lapply(farg, FUN = rep, times = la_m[1L])
    names(arg)[id_fa] <- anm[id_fa]
  }
  y <- do.call(cbind, args = arg) # let warn/err
  if (length(w) != dim(y)[1L] || any(w <= 0)) stop('mixing proportion do not match number of components')
  return(new(Class = 'fmx', parM = y, w = w/sum(w), distname = distname)) # do \strong{not} need to sort by location parameter
}


# parameter names of distributions ('norm' and 'GH' supported)
dist_anm <- function(distname) {
  switch(distname, norm = c('mean', 'sd'), GH = c('A', 'B', 'g', 'h'), stop(sQuote(distname), ' not supported yet'))
}



#' @rdname fmx
#' @export
dfmx <- function(x, dist, distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w, ..., log = FALSE) {
  if (K == 1L) { # no mixture required!!
    switch(distname, 
           norm = return(dnorm(x = x, mean = parM[,1L], sd = parM[,2L], log = log)), 
           GH = return(.dGH(x = x, A = parM[,1L], B = parM[,2L], g = parM[,3L], h = parM[,4L], log = log, ...)),
           stop('distribution ', sQuote(distname), ' not ready yet'))
  }
  
  xm <- tcrossprod(rep(1, times = K), x)
  lds <- switch(distname, 
                norm = dnorm(x = xm, mean = parM[,1L], sd = parM[,2L], log = TRUE), 
                GH = .dGH(x = xm, A = parM[,1L], B = parM[,2L], g = parM[,3L], h = parM[,4L], log = TRUE, ...),
                stop('distribution ', sQuote(distname), ' not ready yet'))
  d <- c(crossprod(w, exp(lds))) # take-log \strong{after} this
  if (!log) return(d)
  
  out <- log(d)
  if (any(id_inf <- is.infinite(out))) {
    stop('this is extremely unlikely to happen')
    # example:
    # (e2 = fmx('GH', A = c(0,3), g = c(.2, .3), h = c(.2, .1), w = c(2, 3)))
    # dfmx(x = -1e20, dist = e2, log = TRUE) # still good
    
    # mixture density being 0, i.e. all component are \approx 0. Use the largest component to replace
    # ds[, id_inf, drop = FALSE] # all 0
    # lds <- lds[, id_inf, drop = FALSE]
    # max_comp <- .Internal(max.col(.Internal(t.default(lds)), 2L))
    # out[id_inf] <- lds[cbind(max_comp, seq_len(sum(id_inf)))] + log(w[max_comp])
  }
  return(out)
}





# not compute intensive
#' @rdname fmx
#' @export
pfmx <- function(q, dist, distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w, ..., lower.tail = TRUE, log.p = FALSE) { # not compute-intensive
  if (K == 1L) {
    switch(distname, 
           norm = return(pnorm(q = q, mean = parM[,1L], sd = parM[,2L], lower.tail = lower.tail, log.p = log.p)),
           GH = return(pGH(q = q, A = parM[,1L], B = parM[,2L], g = parM[,3L], h = parM[,4L], lower.tail = lower.tail, log.p = log.p, ...)),
           stop('distribution ', sQuote(distname), ' not ready yet'))
  }
  
  switch(distname, norm = {
    zM <- tcrossprod(1/parM[,2L], q) - parM[,1L]/parM[,2L]
  }, GH = {
    zM <- tcrossprod(1/parM[,2L], q) - parM[,1L]/parM[,2L]
    gs <- parM[,3L]
    hs <- parM[,4L]
    for (i in seq_len(K)) zM[i,] <- qGH2z(q0 = zM[i,], g = gs[i], h = hs[i], ...)
  }, stop('distribution ', sQuote(distname), ' not ready yet'))
  
  ps <- pnorm(q = zM, lower.tail = lower.tail)
  p <- c(crossprod(w, ps))
  if (log.p) return(log(p)) else return(p)
}


# obtain the `interval` for ?rstpm2::vuniroot
# @rdname fmx
# @export
qfmx_interval <- function(dist, p = c(1e-6, 1-1e-6), distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w, ...) {
  qfun <- paste0('q', distname)
  y_ls <- lapply(seq_len(K), FUN = \(i) {# single component
    iw <- w[i]
    ip <- p
    ip[1L] <- min(p[1L]/iw, .05)
    ip[2L] <- max(1-(1-p[2L])/iw, .95)
    #ip[1L] <- max(p[1L]/iw, .001)
    #ip[2L] <- min(1-(1-p[2L])/iw, .999)
    out <- do.call(qfun, args = c(list(p = ip), as.list.default(parM[i,])))
    if (any(is.infinite(out))) return(invisible()) # very likely to be a mal-estimate
    return(out)
  })
  y <- unlist(y_ls, use.names = FALSE)
  if (anyNA(y)) stop(sQuote(qfun), ' returns NA_real_')
  c(min(y), max(y))
}



#' @rdname fmx
#' @export
qfmx <- function(p, dist, distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w, interval = qfmx_interval(dist = dist), ..., lower.tail = TRUE, log.p = FALSE) {
  # if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval)) stop('masked to save time')
  if (log.p) p <- exp(p)
  
  if (K == 1L) {
    switch(distname, 
           norm = return(qnorm(p, mean = parM[,1L], sd = parM[,2L], lower.tail = lower.tail, log.p = FALSE)),
           GH = return(qGH(p, A = parM[,1L], B = parM[,2L], g = parM[,3L], h = parM[,4L], lower.tail = lower.tail, log.p = FALSE)),
           stop('distribution ', sQuote(distname), ' not ready yet'))
  }
  
  t_w <- t.default(w) # so that not to use `c(.Internal(crossprod(w, ps)))` for CRAN
  #interval <- t.default(array(interval, dim = c(2L, length(p))))
  switch(distname, norm =, GH = {
    sdinv <- 1 / parM[,2L] # constant, save time in zeroin2 algorithm
    eff <- parM[,1L] * sdinv # effect size
  })
  
  switch(distname, norm = {
    #return(vuniroot(f = \(q) { # essentially \code{\link{pfmx}}
    #  zM <- tcrossprod(sdinv, q) - eff # remove .Internal for CRAN
    #  ret <- t_w %*% pnorm(q = zM, lower.tail = lower.tail)
    #  ret - p
    #}, interval = interval)$root)
    return(vuniroot2(y = p, f = \(q) { # essentially \code{\link{pfmx}}
      zM <- tcrossprod(sdinv, q) - eff # remove .Internal for CRAN
      ret <- t_w %*% pnorm(q = zM, lower.tail = lower.tail)
    }, interval = interval))
  }, GH = {
    gs <- parM[,3L]
    hs <- parM[,4L]
    seqid <- seq_len(K)
    #return(vuniroot(f = \(q) { # essentially \code{\link{pfmx}}
    #  z <- q0 <- tcrossprod(sdinv, q) - eff # remove .Internal for CRAN
    #  for (i in seqid) z[i,] <- qGH2z(q0 = q0[i,], g = gs[i], h = hs[i], ...)
    #  ret <- t_w %*% pnorm(q = z, lower.tail = lower.tail)
    #  ret - p 
    #}, interval = interval)$root)
    return(vuniroot2(y = p, f = \(q) { # essentially \code{\link{pfmx}}
      z <- q0 <- tcrossprod(sdinv, q) - eff # remove .Internal for CRAN
      for (i in seqid) z[i,] <- qGH2z(q0 = q0[i,], g = gs[i], h = hs[i], ...)
      ret <- t_w %*% pnorm(q = z, lower.tail = lower.tail)
    }, interval = interval))
    
  }, stop('distribution ', sQuote(distname), ' not ready yet'))
  
}






#' @rdname fmx
#' @export
rfmx <- function(n, dist, distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w) {
  if (!is.integer(n) || anyNA(n) || length(n) != 1L || n <= 0L) stop('sample size must be len-1 positive integer')
  id <- sample.int(n = K, size = n, replace = TRUE, prob = w)
  d2 <- cbind(parM, n = tabulate(id)) # 'matrix'
  r_fn <- paste0('r', distname)
  xs <- lapply(seq_len(K), FUN = \(i) {
    do.call(what = r_fn, args = as.list.default(d2[i, ]))
  })
  out <- unlist(xs, use.names = FALSE)
  return(out)
}


transLog <- function(distname) {
  switch(distname, norm = {
    2L # `sdlog` (no constraint)
  }, GH = {
    c(2L, 4L) # `Blog` & `hlog` (no constraint)
  }, stop('distribution', sQuote(distname), 'not supported yet'))
}





# ?base::`[`
#' @export
`[.fmx` <- function(x, i) {
  y <- x@parM[i, , drop = FALSE]
  w <- x@w
  w <- w / sum(w) # must adjust mixing proportions!!!
  new(Class = 'fmx', parM = y, w = w, distname = x@distname)
}




# extending ?base::as.double
# \strong{not} compute-intensive!
# convert initial values to `par` argument of ?stats::optim
#' @export
as.double.fmx <- function(x, Kmax = K, value_in_name = FALSE, ...) { 
  distname <- x@distname
  parM <- x@parM
  K <- dim(parM)[1L]
  w <- x@w
  if (!is.integer(Kmax) || length(Kmax) != 1L || anyNA(Kmax) || (Kmax < K)) stop('Kmax must be >=', K)
  
  # colnames will be shown as LaTeX in ?knitr::kable
  
  value_in_name <- value_in_name & (Kmax == K)
  
  # plain proportion ('w1' will be removed)
  w_val <- rep(NA_real_, times = Kmax - 1L) # Kmax = 1L compatible
  if (Kmax > 1L) names(w_val) <- if (value_in_name) sprintf('$w_%d=%.0f%%$', 2:Kmax, 1e2*w[2:K]) else sprintf(fmt = '$w_%d$', 2:Kmax)
  if (K > 1L) w_val[1:(K-1L)] <- w[2:K]
   
  # non-log-transformed arguments
  npar <- dim(parM)[2L]
  argnm <- dist_anm(distname)
  if (length(argnm) != npar) stop('typo?')
  dM_new <- if (Kmax == K) parM else rbind(parM, array(NA_real_, dim = c(Kmax-K, npar))) # remove .Internal for publishing on CRAN; not computate-intensive
  out <- c(dM_new, w_val)
  if (value_in_name) {
    names(out)[1:(K*npar)] <- sprintf('$%s_%d=%.1f$', rep(argnm, each = K), seq_len(K), out[1:(K*npar)])
  } else names(out)[1:(Kmax*npar)] <- sprintf('$%s_%d$', rep(argnm, each = Kmax), seq_len(Kmax))
  return(out)
}







# ?base::print
#' @export
print.fmx <- function(x, ...) {
  parM <- x@parM
  K <- dim(parM)[1L]
  parM[] <- sprintf(fmt = '%.2f', parM)
  obj <- if (K == 1L) parM else cbind(parM, w = sprintf(fmt = '%.1f%%', x@w*1e2))
  heading <- paste0(K, '-Component Mixture of ', x@distname, ' Distribution')
  dimnames(obj)[[1L]] <- paste0(seq_len(K), '-comp.')
  cat('\n ', heading, '\n\n', sep = '')
  print.default(obj, quote = FALSE)
  cat('\n')
  print(autoplot.fmx(x))
  cat('\n')
  return(invisible(x))
}


#' @rdname show_S4
#' @export
setMethod(show, signature(object = 'fmx'), definition = \(object) {
  print(object) # \code{print.fmx_QLMDe} or \code{print.fmx}
  return(invisible())
})









