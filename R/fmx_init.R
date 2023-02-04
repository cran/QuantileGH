

#' @title Clustering Observations for Creation of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, one-dimensional observations
#' 
#' @param K \link[base]{integer} scalar, number of mixture components
#' 
#' @param method \link[base]{character} scalar, only \code{'reassign_tkmeans'} supported yet
#' 
#' @param alpha \link[base]{numeric} scalar, proportion of observations to be trimmed in 
#' trimmed k-means algorithm \link[tclust]{tkmeans}
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @return
#' \link{clusterList} returns a \link[base]{list} of \link[base]{numeric} \link[base]{vector}s.
#' 
#' @seealso \link[stats]{kmeans} \link[tclust]{tkmeans} \link{reAssign}
#' 
#' @export
clusterList <- function(x, K, method = c('reassign_tkmeans'), alpha = .05, ...) {
  if (anyNA(x)) stop('input must be free of missing data')
  if (K == 1L) return(list(x))
  tkm <- reAssign.tkmeans(tkmeans(x, k = K, alpha = alpha))
  if (any(1L == tkm$size)) stop('single-observation cluster should be avoided by tclust::tkmeans')
  clus <- tkm[['cluster']]
  attr(clus, which = 'levels') <- as.character.default(seq_len(tkm$k))
  class(clus) <- 'factor'
  ret <- split.default(x, f = clus)[order(c(tkm$centers))] # order by center, only available in dim-1 data!
  return(ret) 
}



#' @title Naive Estimates of Finite Mixture Distribution via Clustering
#' 
#' @description 
#' 
#' Naive estimates for finite mixture distribution \linkS4class{fmx} via clustering.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, observations
#' 
#' @param K \link[base]{integer} scalar, number of mixture components
#' 
#' @param distname \link[base]{character} scalar, name of parametric distribution of the mixture components
#' 
#' @param constraint \link[base]{character} \link[base]{vector}, 
#' parameters (\eqn{g} and/or \eqn{h} for Tukey's \eqn{g}-&-\eqn{h} mixture) to be set at 0.  
#' See \link{fmx_constraint} for details.
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' First of all, if the specified number of components \eqn{K\geq 2},
#' trimmed \eqn{k}-means clustering with re-assignment will be performed;
#' otherwise, all observations will be considered as one single cluster.
#' The standard \eqn{k}-means clustering is not used since the heavy tails of 
#' Tukey's \eqn{g}-&-\eqn{h} distribution could be mistakenly classified as individual cluster(s).
#' 
#' In each of the one or more clusters,
#' \itemize{
#' \item{The letter-value based estimates of Tukey's \eqn{g}-&-\eqn{h} distribution (Hoaglin, 2006)
#' are calculated, for any \eqn{K\geq 1}, serving as the starting values for 
#' QLMD algorithm.   
#' These estimates are provided by \link{fmx_cluster}.}
#' 
#' \item{the \link[stats]{median} and \link[stats:mad]{MAD} will serve as 
#' the starting values for \eqn{\mu} and \eqn{\sigma} 
#' (or \eqn{A} and \eqn{B} for Tukey's \eqn{g}-&-\eqn{h} distribution, with \eqn{g = h = 0}),
#' for QLMD algorithm
#' when \eqn{K = 1}.}
#' }
#' 
#' @return 
#' 
#' \link{fmx_cluster} returns an \linkS4class{fmx} object.
#' 
#' @seealso \link{letterValue}
#' 
#' @export
fmx_cluster <- function(
    x, K, 
    distname = c('GH', 'norm', 'sn'), constraint = character(), 
    ...
) {
  
  xs <- clusterList(x = x, K = K, ...)
  
  dargs <- switch(distname <- match.arg(distname), norm = {
    lapply(xs, FUN = function(i) {
      c(mean = median.default(i), sd = mad(i))
    })
  }, GH = {
    id_constr <- user_constraint(constraint, distname = 'GH', K = K)
    gid <- attr(id_constr, which = 'gid', exact = TRUE)
    hid <- attr(id_constr, which = 'hid', exact = TRUE)
    if (!length(gid) && !length(hid)) {
      if (K == 2L) {
        list(letterValue(xs[[1L]], halfSpread = 'lower'),
             letterValue(xs[[2L]], halfSpread = 'upper'))
      } else lapply(xs, FUN = letterValue)
    } else {
      lapply(seq_len(K), FUN = function(i) {
        ag <- list(p_g = if (i %in% gid) FALSE, p_h = if (i %in% hid) FALSE)
        do.call(letterValue, args = c(list(x = xs[[i]]), ag[lengths(ag, use.names = FALSE) > 0L]))
      })
    }
  }, sn = {
    stop('not programmed')
  }, stop('unknown distname ', distname))

  new(Class = 'fmx', 
      pars = do.call(rbind, args = dargs), 
      w = lengths(xs, use.names = FALSE) / length(x), 
      distname = distname)
}




#' @title Naive Estimates of Finite Mixture Distribution using Mixture of Normal
#' 
#' @description 
#' 
#' Naive estimates for finite mixture distribution \linkS4class{fmx} using mixture of normal
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, observations
#' 
#' @param K \link[base]{integer} scalar, number of mixture components
#' 
#' @param distname \link[base]{character} scalar, name of parametric distribution of the mixture components
#' 
#' @param alpha \link[base]{numeric} scalar, proportion of observations to be trimmed in 
#' trimmed k-means algorithm \link[tclust]{tkmeans}
#' 
#' @param R \link[base]{integer} scalar, number of \link[mixtools]{normalmixEM} replicates
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details
#' 
#' \link{fmx_normix} ... the cluster centers are provided as the starting values of \eqn{\mu}'s for 
#' the univariate normal mixture by EM \link[mixtools:normalmixEM]{algorithm}.
#' \code{R} replicates of normal mixture estimates are obtained, and 
#' the one with maximum likelihood will be selected
#' 
#' @return 
#' 
#' \link{fmx_normix} returns an \linkS4class{fmx} object.
#' 
#' @seealso \link[tclust]{tkmeans} \link[mixtools]{normalmixEM}
#' 
#' @importFrom tclust tkmeans
#' @importFrom mixtools normalmixEM
#' @importFrom utils capture.output
#' @importFrom stats median.default
#' 
#' @export
fmx_normix <- function(x, K, distname = c('GH', 'norm', 'sn'), alpha = .05, R = 10L, ...) {
  
  if (anyNA(x)) stop('input must be free of missing data')
  distname <- match.arg(distname)
  
  if (K == 1L) {
    mu <- median.default(x)
    sigma <- mad(x)
    return(new(Class = 'fmx', distname = distname, pars = switch(distname, norm = cbind(
      mean = mu, sd = sigma
    ), GH = cbind(
      A = mu, B = sigma, g = 0, h = 0
    ), sn = cbind(
      xi = mu, omega = sigma, alpha = 0  # I dont understand `tau`
    ), stop('`distname` not supported ', distname))))
  }
  
  tkm <- tkmeans(x, k = K, alpha = alpha)
  tx <- x[tkm$cluster != 0] # use un-trimmed observations!!!
  
  rets <- replicate(n = R, expr = {
    invisible(capture.output(tmp <- normalmixEM(
      x = tx, k = K, 
      mu = sort.int(c(tkm$centers)), # hope location-pars not get too close..
      lambda = rep(1/K, times = K), # hope none of mix-prop-pars gets too small ..
      epsilon = 1e-2)))
    tmp
  }, simplify = FALSE)
  
  tmp <- sort.mixEM(rets[[which.max(vapply(rets, FUN = logLik.mixEM, FUN.VALUE = 0, USE.NAMES = FALSE))]])
  
  return(new(Class = 'fmx', distname = distname, w = tmp$lambda, pars = switch(distname, norm = cbind(
    mean = tmp$mu, sd = tmp$sigma
  ), GH = cbind(
    A = tmp$mu, B = tmp$sigma, g = 0, h = 0
  ), sn = cbind(
    xi = tmp$mu, omega = tmp$sigma, alpha = 0 # I dont understand `tau`
  ))))
  
}







#' @title Best Naive Estimates for Finite Mixture Distribution
#' 
#' @description 
#' 
#' Best estimates for finite mixture distribution \linkS4class{fmx}.
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, observations
#' 
#' @param test \link[base]{character} scalar, criteria for selecting the optimal estimates.  
#' See \strong{Details}.
#' 
#' @param ... additional parameters of \link{fmx_normix} and \link{fmx_cluster}
#' 
#' @details
#' 
#' \link{fmx_hybrid} compares 
#' the Tukey's \eqn{g}-&-\eqn{h} mixture estimate provided by \link{fmx_cluster}
#' and the normal mixture estimate by \link{fmx_normix}, 
#' and select the one either with maximum likelihood (\code{test = 'logLik'}, default), 
#' with minimum Cramer-von Mises distance (\code{test = 'CvM'}) or 
#' with minimum Kolmogorov distance (\code{test = 'Kolmogorov'}).
#' 
#' @return 
#' 
#' \link{fmx_hybrid} returns an \linkS4class{fmx} object.
#' 
#' @examples 
#' d1 = fmx('norm', mean = c(1, 2), sd = .5, w = c(.4, .6))
#' set.seed(100); hist(x1 <- rfmx(n = 1e3L, dist = d1))
#' fmx_normix(x1, distname = 'norm', K = 2L)
#' fmx_normix(x1, distname = 'GH', K = 2L)
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(100); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' fmx_cluster(x2, K = 2L)
#' fmx_cluster(x2, K = 2L, constraint = c('g1', 'h2'))
#' fmx_normix(x2, K = 2L, distname = 'GH')
#' fmx_hybrid(x2, distname = 'GH', K = 2L)
#' 
#' @export
fmx_hybrid <- function(x, test = c('logLik', 'CvM', 'KS'), ...) {
  
  test <- match.arg(test)
  
  ys <- list(
    cluster = fmx_cluster(x, ...),
    normix = fmx_normix(x, ...)
  )

  stats <- list(
    logLik = - vapply(ys, FUN = logLik.fmx, data = x, FUN.VALUE = NA_real_),
    CvM = vapply(ys, FUN = function(i) fmx_diagnosis(i, data = x, type = 'CramerVonMises')$statistic, FUN.VALUE = NA_real_),
    Kolmogorov = vapply(ys, FUN = function(i) fmx_diagnosis(i, data = x, type = 'Kolmogorov'), FUN.VALUE = NA_real_)
  )
  if (anyNA(stats, recursive = TRUE)) stop('should not happen')
  
  chosen <- vapply(stats, FUN = function(i) names(which.min(i)), FUN.VALUE = '')
  
  ret <- ys[[chosen[test]]]
  if (TRUE) { # for developer
    #attr(ret, 'ys') <- ys
    #attr(ret, 'x') <- x
    #attr(ret, 'fig') <- batchplot_fmx(ys, obs = x) # slow
    #attr(ret, 'chosen') <- chosen
  }
  return(ret)
}




#stop('make it clear in the documentation that the user can choose to use lettervalue or normix')

#stop('make it clear in the documentation that we offer the hybrid of lettervalue and normix')
