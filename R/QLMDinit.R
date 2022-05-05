



#' @title Initial Values for Quantile Least Mahalanobis Distance (QLMD) Estimates 
#' 
#' @description 
#' 
#' Various methods for obtaining the initial values for 
#' Quantile Least Mahalanobis Distance (QLMD)
#' estimates of finite mixture distribution \linkS4class{fmx}.
#' 
#' @param x \link[base]{numeric} vector, the actual observations
#' 
#' @param distname \link[base]{character} scalar, name of parametric distribution
#' 
#' @param K \link[base]{integer} scalar, number of mixture components
#' 
#' @param constraint \link[base]{character} vector, 
#' parameters (\eqn{g} and/or \eqn{h} for Tukey's \eqn{g}-&-\eqn{h} mixture) to be set at 0.  
#' See \link{fmx_constraint} for details.
#' 
#' @param alpha \link[base]{numeric} scalar, proportion of observations to be trimmed by 
#' trimmed k-means algorithm \link[tclust]{tkmeans}
#' 
#' @param R \link[base]{integer} scalar, number of \link[mixtools]{normalmixEM} replicates
#' 
#' @param test \link[base]{character} scalar, criteria for selecting the optimal estimates.  
#' See \strong{Details}.
#' 
#' @param ... additional arguments
#' 
#' @details
#' 
#' First of all, if the specified number of components \eqn{K\geq 2},  
#' trimmed \eqn{k}-means clustering with re-assignment will be performed;
#' otherwise, all observations will be considered as one single cluster.
#' The standard \eqn{k}-means clustering is not used since the heavy tails of 
#' Tukey's \eqn{g}-&-\eqn{h} distribution could be mistakenly classified as individual cluster(s).
#' In each of the one or more clusters,
#' 
#' \itemize{
#' 
#' \item{The letter-value based estimates of Tukey's \eqn{g}-&-\eqn{h} distribution (Hoaglin, 2006)
#' are calculated, for any \eqn{K\geq 1}, serving as the starting values for 
#' QLMD algorithm.   
#' These estimates are provided by \code{\link{QLMDinit_letterValue}}.}
#' 
#' \item{the \link[stats]{median} and \link[stats:mad]{MAD} will serve as 
#' the starting values for \eqn{\mu} and \eqn{\sigma} 
#' (or \eqn{A} and \eqn{B} for Tukey's \eqn{g}-&-\eqn{h} distribution, with \eqn{g = h = 0}),
#' for QLMD algorithm
#' when \eqn{K = 1}.
#' Otherwise, the cluster centers are provided as the starting values of \eqn{\mu}'s for 
#' the univariate normal mixture by EM \link[mixtools:normalmixEM]{algorithm}.
#' \code{R} replicates of normal mixture estimates are obtained, and 
#' the one with maximum likelihood will serve as the starting values for 
#' QLMD algorithm.
#' These estimates are provided by \code{\link{QLMDinit_normix}}.}
#' 
#' }
#' 
#' \code{\link{QLMDinit}} compares 
#' the Tukey's \eqn{g}-&-\eqn{h} mixture estimate provided by \code{\link{QLMDinit_letterValue}}
#' and the normal mixture estimate by \code{\link{QLMDinit_normix}}, 
#' and select the one either with maximum likelihood (\code{test = 'logLik'}, default), 
#' with minimum Cramer-von Mises distance (\code{test = 'CvM'}) or 
#' with minimum Kolmogorov-Smirnov distance (\code{test = 'KS'}).
#' 
#' @return 
#' 
#' \code{\link{QLMDinit_letterValue}}, \link{QLMDinit_normix} and \link{QLMDinit}
#' all return \linkS4class{fmx} objects.
#' 
#' @seealso \link[stats]{kmeans} \link[tclust]{tkmeans} \link{reAssign} \link{letterValue} \link[mixtools]{normalmixEM}
#' 
#' @examples 
#' d1 = fmx('norm', mean = c(1, 2), sd = .5, w = c(.4, .6))
#' set.seed(100); hist(x1 <- rfmx(n = 1e3L, dist = d1))
#' QLMDinit_normix(x1, distname = 'norm', K = 2L)
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(100); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' QLMDinit_letterValue(x2, K = 2L)
#' QLMDinit_letterValue(x2, K = 2L, constraint = c('g1', 'h2'))
#' QLMDinit_normix(x2, K = 2L)
#' QLMDinit(x2, distname = 'GH', K = 2L)
#' 
#' @name QLMDinit
#' @export
QLMDinit_letterValue <- function(x, K, constraint = character(), alpha = .05, ...) {
  if (K == 1L) {
    w <- 1
    xs <- list(x)
  } else {
    # either ?stats::kmeans or ?tclust::tkmeans are slow for big `x`
    tkm <- reAssign.tkmeans(tkmeans(x, k = K, alpha = alpha))
    if (any(1L == tkm$size)) stop('single-observation cluster should be avoided by tclust::tkmeans')
    clus <- tkm[['cluster']]
    attr(clus, which = 'levels') <- as.character.default(seq_len(tkm$k))
    class(clus) <- 'factor'
    # order by center, only available in dim-1 data!
    o <- order(c(tkm$centers)) # not compute-intensive
    xs0 <- split.default(x, f = clus) # identical to .Internal(split(x, clus)); not compute intensive
    xs <- xs0[o]
    w <- tkm$weights[o]
  }
  
  id_constr <- fmx_constraint_user(distname = 'GH', K = K, user = constraint)
  gid <- attr(id_constr, which = 'gid', exact = TRUE)
  hid <- attr(id_constr, which = 'hid', exact = TRUE)
  dargs <- if (!length(gid) && !length(hid)) {
    if (K == 2L) {
      list(letterValue(xs[[1L]], halfSpread = 'lower'),
           letterValue(xs[[2L]], halfSpread = 'upper'))
    } else lapply(xs, FUN = letterValue)
  } else {
    lapply(seq_len(K), FUN = function(i) {
      ag <- list(
        p_g = if (i %in% gid) FALSE,
        p_h = if (i %in% hid) FALSE
      )
      do.call(letterValue, args = c(list(x = xs[[i]]), ag[lengths(ag, use.names = FALSE) > 0L]))
    })
  }

  return(new(Class = 'fmx', parM = do.call(rbind, args = dargs), w = w, distname = 'GH'))
}


#' @rdname QLMDinit
#' @export
QLMDinit_normix <- function(x, K, alpha = .05, R = 10L, ...) {
  
  if (K == 1L) return(new(Class = 'fmx', distname = 'norm', parM = cbind(mean = median.default(x), sd = mad(x)))) 
  
  tkm <- tkmeans(x, k = K, alpha = alpha)
  tx <- x[tkm$cluster != 0]
  
  rets <- replicate(n = R, expr = {
    invisible(capture.output(tmp <- normalmixEM(
      x = tx, k = K, 
      mu = sort.int(c(tkm$centers)), # hope location-pars not get too close..
      lambda = rep(1/K, times = K), # hope none of mix-prop-pars gets too small ..
      epsilon = 1e-2)))
    tmp
  }, simplify = FALSE)
  
  tmp <- sort.mixEM(rets[[which.max(vapply(rets, FUN = logLik.mixEM, FUN.VALUE = 0, USE.NAMES = FALSE))]])
  
  ret <- new(Class = 'fmx', parM = cbind(mean = tmp$mu, sd = tmp$sigma), w = tmp$lambda, distname = 'norm')
  #, GH = cbind(
  #  A = tmp$mu, B = tmp$sigma, g = 0, h = 0
  #)), 
  
  #attr(ret, which = 'mixEM') <- tmp
  return(ret)
}


#' @rdname QLMDinit
#' @export
QLMDinit <- function(x, distname = c('GH', 'norm'), test = c('logLik', 'CvM', 'KS'), ...) {
  
  distname <- match.arg(distname)
  test <- match.arg(test)
  
  y_norm <- QLMDinit_normix(x, ...)
  if (distname == 'norm') return(y_norm)
  
  inits <- list(
    letterValue = QLMDinit_letterValue(x, ...),
    normix = new(Class = 'fmx', parM = cbind(y_norm@parM, 0, 0), w = y_norm@w, distname = 'GH')
  )

  stats <- list(
    logLik = - vapply(inits, FUN = logLik.fmx, data = x, FUN.VALUE = NA_real_),
    CvM = vapply(inits, FUN = function(i) CvM_test.fmx(i, data = x)$statistic, FUN.VALUE = NA_real_),
    KS = vapply(inits, FUN = function(i) ks_test.fmx(i, data = x)$statistic, FUN.VALUE = NA_real_)
  )
  if (anyNA(stats, recursive = TRUE)) stop('do not allow NA in resuls of logLik, CvM_test and ks_test')
  
  chosen <- vapply(stats, FUN = function(i) names(which.min(i)), FUN.VALUE = '')
  
  ret <- inits[[chosen[test]]]
  if (TRUE) { # for developer
    #attr(ret, 'inits') <- inits
    #attr(ret, 'x') <- x
    #attr(ret, 'fig') <- batchplot_fmx(inits, obs = x) # slow
    attr(ret, 'chosen') <- chosen
  }
  return(ret)
}



# see ?mixtools::normalmixEM carefully
# ?stats::logLik
#' @export
logLik.mixEM <- function(object, ...) {
  val <- object[['loglik']]
  attr(val, which = 'nobs') <- length(object[['x']])
  attr(val, which = 'df') <- switch(object[['ft']], normalmixEM = {
    3 * length(object[['mu']]) - 1L
  }, stop('what?'))
  class(val) <- 'logLik'
  return(val)
}

# ?base::sort
#' @export
sort.mixEM <- function(x, decreasing = FALSE, ...) {
  o <- order(x[['mu']], decreasing = decreasing)
  # ?mixtools::normalmixEM does \strong{not} order the location parameter
  y <- x
  y[['lambda']] <- x[['lambda']][o]
  y[['mu']] <- x[['mu']][o]
  y[['sigma']] <- x[['sigma']][o]
  y[['posterior']] <- x[['posterior']][, o]
  colnames(y[['posterior']]) <- paste0('comp.', seq_along(x[['mu']]))
  return(y)
}


