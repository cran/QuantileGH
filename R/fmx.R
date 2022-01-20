

#' @title The Finite Mixture Distribution
#' 
#' @description 
#' 
#' Density function, distribution function, quantile function and random generation for a finite mixture distribution 
#' with normal or Tukey's \eqn{g}-&-\eqn{h} components.
#' 
#' @param x,q vector of quantiles, \code{NA_real_} value(s) allowed.
#' 
#' @param p vector of probabilities.
#' 
#' @param n number of observations.
#' 
#' @param dist \code{\linkS4class{fmx}} object, representing a finite mixture distribution
#' 
#' @param silent \code{\link[base]{logical}} scalar, whether to print out necessary messages, default \code{TRUE}
#' 
#' @param log,log.p \code{\link[base]{logical}} scalar; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#' 
#' @param lower.tail \code{\link[base]{logical}} scalar; if \code{TRUE} (default), probabilities are \eqn{Pr(X\le x)} otherwise, \eqn{Pr(X>x)}.
#' 
#' @param interval interval for root finding (see \code{\link[rstpm2]{vuniroot}})
#' 
#' @param distname,K,parM,w auxiliary parameters, whose default values are determined by
#' the \code{\linkS4class{fmx}} object provided in argument \code{dist}.
#' The user-specified vector of \code{w} does not need to sum up to 1; \code{w/sum(w)} will be used internally.
#' 
#' @param ... mixture distribution parameters in \code{\link{fmx}}.
#' See \code{\link{dGH}} for the names and default values of Tukey's \eqn{g}-&-\eqn{h} distribution parameters, 
#' or \code{\link[stats]{dnorm}} for the names and default values of normal distribution parameters.
#' 
#' @details 
#' 
#' A computational challenge in \code{\link{dfmx}} is when mixture density is very close to 0,
#' which happens when the per-component log densities are negative with big absolute values.  
#' In such case, we cannot compute the log mixture densities (i.e., \code{-Inf}), for the log-likelihood \code{\link{logLik.fmx_QLMDe}}.
#' Our solution is to replace these \code{-Inf} log mixture densities by 
#' the weighted average (using the mixing proportions of \code{dist}) 
#' of the per-component log densities.
#' 
#' \code{\link{qfmx}} gives the quantile function, by solving \code{\link{pfmx}} by \code{\link[rstpm2]{vuniroot}}.
#' One major challenge when dealing with the finite mixture of Tukey's \eqn{g}-&-\eqn{h} family distribution
#' is that Brent–Dekker's method needs to be performed in both \code{\link{pGH}} and \code{\link{qfmx}}, 
#' i.e. `two layers' of root-finding algorithm.
#' 
#' 
#' @return 
#' 
#' \code{\link{fmx}} returns an \code{\linkS4class{fmx}} object which specifies the parameters of a finite mixture distribution.
#' 
#' \code{\link{dfmx}} returns a vector of probability density values of an \code{\linkS4class{fmx}} object at specified quantiles \code{x}.
#' 
#' \code{\link{pfmx}} returns a vector of cumulative probability values of an \code{\linkS4class{fmx}} object at specified quantiles \code{q}.
#' 
#' \code{\link{qfmx}} returns an unnamed vector of quantiles of an \code{\linkS4class{fmx}} object, based on specified cumulative probabilities \code{p}.
#' Note that \code{\link[stats]{qnorm}} returns an unnamed vector of quantiles, 
#' although \code{\link[stats]{quantile}} returns a named vector of quantiles.
#' 
#' \code{\link{rfmx}} generates random deviates of an \code{\linkS4class{fmx}} object.
#' 
#' @examples 
#' 
#' # paramter is recycled
#' fmx('norm', mean = c(4, 1, 14, 11), w = c(1, 2))
#' 
#' x = (-3):7
#' 
#' (e1 = fmx('norm', mean = c(0,3), sd = c(1,1.3), w = c(1, 1)))
#' isS4(e1) # TRUE
#' slotNames(e1)
#' plot(e1) # using vanilla R
#' autoplot(e1) # using ggplot2 package
#' hist(rfmx(n = 1e3L, dist = e1), main = '1000 obs from e1')
#' # generate a sample of size 1e3L from mixture distribution `e1`
#' round(dfmx(x, dist = e1), digits = 3L)
#' round(p1 <- pfmx(x, dist = e1), digits = 3L)
#' stopifnot(all.equal.numeric(qfmx(p1, dist = e1), x, tol = 1e-4))
#' 
#' (e2 = fmx('GH', A = c(0,3), g = c(.2, .3), h = c(.2, .1), w = c(2, 3)))
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
#' if (FALSE) {
#'   # log-mixture-density smoothing, for developers
#'   (e4 = fmx('norm', mean = c(0,3), w = c(2, 3)))
#'   plot(e4, type = 'log_pdf', xlim = c(-50, 50))
#' }
#' 
#' @name fmx
#' @export
fmx <- function(distname, w = 1, ...) {
  if (!is.character(distname) || (length(distname) != 1L) || anyNA(distname) || !nzchar(distname)) stop('distname must be len-1 char')
  if (any(w <= 0)) stop('mixing proportion do not match number of components')
  
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
  
  if (anyNA(arg, recursive = TRUE)) stop('do not allow NA in `arg`')
  y <- do.call(cbind, args = c(arg, list(w))) # let warn/err (long is not a multiple of short, etc.)
  if (is.unsorted(loc <- y[,1L], strictly = FALSE)) {
    message('Re-ordered by location parameter\n')
    y <- y[order(loc), , drop = FALSE]
  }
  dy <- dim(y)
  n <- dy[2L]
  new(Class = 'fmx', parM = y[, -n, drop = FALSE], w = unname(y[,n]/sum(y[,n])), distname = distname)
}


# parameter names of distributions ('norm' and 'GH' supported)
dist_anm <- function(distname) {
  switch(distname, norm = c('mean', 'sd'), GH = c('A', 'B', 'g', 'h'), stop(sQuote(distname), ' not supported yet'))
}



#' @rdname fmx
#' @export
dfmx <- function(x, dist, distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w, silent = TRUE, ..., log = FALSE) {
  if (K == 1L) { # no mixture required!!
    switch(distname, 
           norm = return(dnorm(x = x, mean = parM[,1L], sd = parM[,2L], log = log)), 
           GH = return(.dGH(x = x, A = parM[,1L], B = parM[,2L], g = parM[,3L], h = parM[,4L], log = log, ...)),
           stop('distribution ', sQuote(distname), ' not ready yet'))
  }
  
  xm <- tcrossprod(rep(1, times = K), x)
  
  lds <- switch(distname, # `lds` is per-component log-densities
                norm = dnorm(x = xm, mean = parM[,1L], sd = parM[,2L], log = TRUE), 
                GH = .dGH(x = xm, A = parM[,1L], B = parM[,2L], g = parM[,3L], h = parM[,4L], log = TRUE, ...),
                stop('distribution ', sQuote(distname), ' not ready yet'))
  if (any(is.infinite(lds))) {
    #tmp <<- dist
    stop('per-component log-density should not be -Inf (unless `sd` or `B` is 0)') # is.infinite(NA) is \code{FALSE}
  }
  
  d <- c(crossprod(w, exp(lds)))
  if (!log) return(d)
  #any(is.infinite(log(d))) # could happen
  
  nx <- length(x)
  wlds <- log(w) + lds # weighted-lds
  id_max <- max.col(t.default(wlds))  # dont want to import ?matrixStats::colMaxs
  id_max_seq <- cbind(id_max, seq_len(nx))
  logd <- wlds[id_max_seq] + log(.colSums(tcrossprod(w, 1/w[id_max]) * exp(t.default(t.default(lds) - lds[id_max_seq])), m = K, n = nx))
  if (!isTRUE(all.equal.numeric(exp(logd), d))) stop('new dfmx wrong?')
  return(logd)
  
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
  y_ls <- lapply(seq_len(K), FUN = function(i) {# single component
    iw <- w[i]
    ip <- p
    ip[1L] <- min(p[1L]/iw, .05)
    ip[2L] <- max(1-(1-p[2L])/iw, .95)
    #ip[1L] <- max(p[1L]/iw, .001)
    #ip[2L] <- min(1-(1-p[2L])/iw, .999)
    out <- do.call(qfun, args = c(list(p = ip), as.list.default(parM[i,])))
    if (any(is.infinite(out))) return(invisible()) # very likely to be a malformed estimate
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
  
  t_w <- t.default(w)
  switch(distname, norm =, GH = {
    sdinv <- 1 / parM[,2L] # constant, save time in vuniroot algorithm
    eff <- parM[,1L] * sdinv # effect size
  })
  
  switch(distname, norm = {
    return(vuniroot2(y = p, f = function(q) { # essentially \code{\link{pfmx}}
      zM <- tcrossprod(sdinv, q) - eff
      c(t_w %*% pnorm(q = zM, lower.tail = lower.tail))
    }, interval = interval))
  }, GH = {
    gs <- parM[,3L]
    hs <- parM[,4L]
    seqid <- seq_len(K)
    return(vuniroot2(y = p, f = function(q) { # essentially \code{\link{pfmx}}
      z <- q0 <- tcrossprod(sdinv, q) - eff
      for (i in seqid) z[i,] <- qGH2z(q0 = q0[i,], g = gs[i], h = hs[i], ...)
      c(t_w %*% pnorm(q = z, lower.tail = lower.tail))
    }, interval = interval))
    
  }, stop('distribution ', sQuote(distname), ' not ready yet'))
  
}






#' @rdname fmx
#' @export
rfmx <- function(n, dist, distname = dist@distname, K = dim(parM)[1L], parM = dist@parM, w = dist@w) {
  if (!is.integer(n) || anyNA(n) || length(n) != 1L || n <= 0L) stop('sample size must be len-1 positive integer (e.g., use 100L instead of 100)')
  id <- sample.int(n = K, size = n, replace = TRUE, prob = w)
  d2 <- cbind(parM, n = tabulate(id)) # 'matrix'
  r_fn <- paste0('r', distname)
  xs <- lapply(seq_len(K), FUN = function(i) {
    do.call(what = r_fn, args = as.list.default(d2[i, ]))
  })
  out <- unlist(xs, use.names = FALSE)
  return(out)
}






#' @title Number of Components in \code{\linkS4class{fmx}} and \code{\linkS4class{fmx_QLMDe}} Object
#' 
#' @description
#' 
#' Obtain the number of components in \code{\linkS4class{fmx}} and \code{\linkS4class{fmx_QLMDe}} object.
#' 
#' @param x \code{\linkS4class{fmx}} and \code{\linkS4class{fmx_QLMDe}} object.
#' 
#' @details 
#' 
#' For user convenience
#' 
#' @return 
#' 
#' An \code{\link[base]{integer}} value indicating the number of components in 
#' an \code{\linkS4class{fmx}} and.or \code{\linkS4class{fmx_QLMDe}} object.
#' 
#' @examples 
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' K.fmx(d2)
#' 
#' @export
K.fmx <- function(x) dim(x@parM)[1L] # only used in high level user-interface




#' @title Subset of Components in \code{\linkS4class{fmx}} and/or \code{\linkS4class{fmx_QLMDe}} Object
#' 
#' @description 
#' 
#' Taking subset of components in \code{\linkS4class{fmx}} and/or \code{\linkS4class{fmx_QLMDe}} object
#' 
#' @param x \code{\linkS4class{fmx}} and/or \code{\linkS4class{fmx_QLMDe}} object
#' 
#' @param i \code{\link[base]{integer}} or \code{\link[base]{logical}} vector, 
#' the row index(es) of the subset of components to be chosen, see \code{\link[base]{[}}
#' 
#' @param j ignored (always \code{TRUE}, i.e., all parameters of such give distribution must be selected), see \code{\link[base]{[}}
#' 
#' @param drop ignored (always \code{FALSE}), see \code{\link[base]{[}}
#' 
#' @details 
#' 
#' Note that using definitions as S3 method dispatch \code{`[.fmx`} or \code{`[.fmx_QLMDe`} won't work 
#' for S4 objects \code{\linkS4class{fmx}} and/or \code{\linkS4class{fmx_QLMDe}}.
#' 
#' @return 
#' 
#' An \code{\linkS4class{fmx}} object consisting of a subset of components.
#' Note that subsetting \code{\linkS4class{fmx_QLMDe}} object will return an \code{\linkS4class{fmx}} object, 
#' which contains only the mixture parameters, i.e., information about the observations (e.g. slots \code{@@data} and \code{@@data.name}),
#' as well as other estimation related slots (e.g., \code{@@init}) will be lost.
#' 
#' @examples 
#' 
#' (d = fmx('norm', mean = c(1, 5, 9)))
#' d[1:2, ]
#' 
#' @export
setMethod(`[`, signature(x = 'fmx', i = 'ANY', j = 'ANY', drop = 'ANY'), definition = function(x, i, j, drop) {
  if (inherits(x, what = 'fmx_QLMDe')) message('Subsetting `fmx_QLMDe` will subset the estimates and drop `@data` etc.')
  parM <- x@parM[i, , drop = FALSE]
  w <- x@w[i]
  w <- unname(w / sum(w)) # adjust mixing proportions
  o <- order(parM[, 1L])
  new(Class = 'fmx', parM = parM[o, , drop = FALSE], w = w[o], distname = x@distname)
})






# ?base::print
#' @export
print.fmx <- function(x, ...) {
  parM <- x@parM
  K <- dim(parM)[1L]
  parM[] <- sprintf(fmt = '%.2f', parM)
  if (length(id_constr <- fmx_constraint(x))) {
    parM[id_constr] <- '.'
  }
  obj <- if (K == 1L) parM else cbind(parM, w = sprintf(fmt = '%.1f%%', x@w*1e2))
  heading <- paste0(K, '-Component Mixture of ', switch(x@distname, norm = 'Normal', GH = 'Tukey\'s G-&-H'), ' Distribution')
  dimnames(obj)[[1L]] <- paste0(seq_len(K), '-comp.')
  cat('\n ', heading, '\n\n', sep = '')
  print.default(obj, quote = FALSE)
  if (length(id_constr)) cat('\nwhere ', sQuote('.'), ' denotes an enforced constraint\n', sep = '')
  cat('\n')
  print(autoplot.fmx(x))
  return(invisible(x))
}


#' @rdname show_S4
#' @export
setMethod(show, signature(object = 'fmx'), definition = function(object) {
  print(object) # \code{\link{print.fmx_QLMDe}} or \code{\link{print.fmx}}
  return(invisible())
})









