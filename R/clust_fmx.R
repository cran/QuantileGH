
#' @title Estimates of Finite Mixture Distribution (\code{linkS4class{fmx}}) via Clustering
#' 
#' @description 
#' 
#' First, trimmed \eqn{k}-means clustering with re-assignment is used to assign all trimmed observations 
#' back to one of the \eqn{K} components, 
#' then for each component, robust parameter estimates are computed for each 
#' mixture component using only observations from the corresponding cluster.
#' 
#' @param x \code{\link[base]{numeric}} vector of observations
#' 
#' @param distname \code{\link[base]{character}} value for the name of parametric distribution
#' 
#' @param K \code{\link[base]{integer}} value of the number of components
#' 
#' @param constraint see \code{\link{QLMDe}}
#' 
#' @return 
#' 
#' An S4 \code{\linkS4class{fmx}} object fitted on given observations, using trimmed \eqn{k}-means clustering with re-assignment,
#' and cluster-wise robust parameter estimates.
#' 
#' @details
#' 
#' The naive estimates for the mixture distribution is obtained in the following two steps.
#' 
#' \itemize{
#' 
#' \item {Trimmed \eqn{k}-means clustering (\code{\link[tclust]{tkmeans}}) with re-assignment. 
#' First, \eqn{5\%} (default value of parameter \code{alpha} of \code{\link[tclust]{tkmeans}}) 
#' of observations are trimmed and \eqn{k}-mean clustering is performed.
#' Next, the \code{\link[stats]{mahalanobis}} distance is computed between each trimmed observation and each cluster.
#' Then, each trimmed observation is assigned to the closest cluster (at the smallest Mahalanobis distance). 
#' The standard \eqn{k}-means clustering (\code{\link[stats]{kmeans}}) is not used since the heavy tails of 
#' Tukey's \eqn{g}-&-\eqn{h} distribution can be mistakenly classified as a cluster.}
#' 
#' \item {In each cluster, the following estimates are computed, 
#' (i). median and median absolute deviation (\code{\link[stats]{mad}}) for individual normal component;
#' (ii). Letter-based estimates for individual Tukey's \eqn{g}-&-\eqn{h} components.}
#' 
#' \item {The mixing probability of each component is estimated by the proportion of observations 
#' in the corresponding cluster.}
#' 
#' }
#' 
#' @references 
#' David C. Hoaglin (2006)  Summarizing Shape Numerically: The \eqn{g}-and-\eqn{h} Distributions. 
#' Wiley Series in Probability and Statistics. Chapter 11.  \doi{10.1002/9781118150702.ch11}
#' 
#' 
#' @examples 
#' d1 = fmx('norm', mean = c(1, 2), sd = .5, w = c(.4, .6))
#' set.seed(100); hist(x1 <- rfmx(n = 1e3L, dist = d1))
#' clust_fmx(x1, distname = 'norm', K = 2L)
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(100); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' clust_fmx(x2, distname = 'GH', K = 2L)
#' clust_fmx(x2, distname = 'GH', K = 2L, constraint = c('g1', 'h2'))
#' 
#' @export
clust_fmx <- function(x, distname, K, constraint = character()) {
  if (K == 1L) {
    w <- 1
    xs <- list(x)
  } else {
    # either ?stats::kmeans or ?tclust::tkmeans are slow for big `x`
    tkm <- reAssign.tkmeans(tkmeans(x, k = K, alpha = .05)) # imported ?tclust::tkmeans
    if (any(1L == tkm$size)) stop('single-observation cluster should be avoided by tclust::tkmeans')
    clus <- tkm[['cluster']]
    attr(clus, 'levels') <- as.character.default(seq_len(tkm$k))
    class(clus) <- 'factor'
    # order by center, only available in dim-1 data!
    o <- order(c(tkm$centers)) # not compute-intensive
    xs0 <- split.default(x, f = clus) # identical to .Internal(split(x, clus)); not compute intensive
    xs <- xs0[o]
    w <- tkm$weights[o]
  }
  
  dargs <- switch(distname, GH = {
    id_constr <- fmx_constraint_user(distname = distname, K = K, user = constraint)
    gid <- attr(id_constr, which = 'gid', exact = TRUE)
    hid <- attr(id_constr, which = 'hid', exact = TRUE)
    if (!length(gid) && !length(hid)) {
      if (K == 2L) {
        list(Hoaglin_GH(xs[[1L]], halfSpread = 'lower'),
             Hoaglin_GH(xs[[2L]], halfSpread = 'upper'))
      } else lapply(xs, FUN = Hoaglin_GH)
    } else {
      lapply(seq_len(K), FUN = function(i) {
        ag <- list(
          p_g = if (i %in% gid) FALSE,
          p_h = if (i %in% hid) FALSE
        )
        do.call(Hoaglin_GH, args = c(list(x = xs[[i]]), ag[lengths(ag, use.names = FALSE) > 0L]))
      })
    }
  }, norm = {
    lapply(xs, FUN = function(ix) {
      x_med <- median.default(ix) # to match Hoaglin's A (not to use \code{c(tkm$centers)[o]} directly)
      return(c(mean = x_med, sd = mad(ix, center = x_med))) # not compute-intensive
    })
  }, stop('distribution not supported yet: ', sQuote(distname)))
  
  return(new(Class = 'fmx', parM = do.call(rbind, args = dargs), w = w, distname = distname))
}



#' @title Parameter Constraint(s) of Mixture Distribution
#' 
#' @description 
#' 
#' Determine the parameter constraint(s) of a finite mixture distribution, either 
#' by the value of parameters of such mixture distribution, or by a user-specified string. 
#' 
#' @param dist an \code{\linkS4class{fmx}} object, can be missing
#' 
#' @param distname,K,parM the name of distribution, the number of components and the matrix of distribution parameters of a finite mixture distribution
#' 
#' @param user an user-specified \code{\link[base]{character}} vector to denote the constraint(s) to be 
#' imposed for a finite mixture distribution.  For example, for a two-component Tukey's \eqn{g}-&-\eqn{h}
#' mixture, \code{user = c('g2', 'h1')} indicates the \eqn{g}-parameter for the first component (with smaller mean value)
#' and the \eqn{h}-parameter for the second component (with larger mean value) are to be constrained, i.e., \eqn{g_2=h_1=0}.
#' 
#' @return 
#' 
#' \code{\link{fmx_constraint}} returns the indexes of internal parameters 
#' (only applicable to Tukey's \eqn{g}-&-\eqn{h} mixture distribution, yet) to be constrained, 
#' based on the input \code{\linkS4class{fmx}} object \code{dist}.
#' 
#' \code{\link{fmx_constraint_user}} returns the indexes of internal parameters 
#' (only applicable to Tukey's \eqn{g}-&-\eqn{h} mixture distribution, yet) to be constrained, 
#' based on the type of distribution (\code{distname}), number of components (\code{K}) 
#' and a user-specified string (e.g., \code{c('g2', 'h1')}).
#' 
#' \code{\link{fmx_constraint_brief}} returns a \code{\link[base]{character}} scalar (of LaTeX expression) of the constraint, 
#' primarily intended for end-users in plots.
#' 
#' 
#' @examples 
#' (d0 = fmx('GH', A = c(1,4), g = c(.2,.1), h = c(.05,.1), w = c(1,1)))
#' (c0 = fmx_constraint(d0))
#' stopifnot(identical(c0, fmx_constraint_user(distname = 'GH', K = 2L, user = character())))
#' fmx_constraint_brief(d0)
#' 
#' (d1 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(0,.1), w = c(1,1)))
#' (c1 = fmx_constraint(d1))
#' stopifnot(identical(c1, fmx_constraint_user(distname = 'GH', K = 2L, user = c('g2', 'h1'))))
#' fmx_constraint_brief(d1)
#' 
#' (d2 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(.15,.1), w = c(1,1)))
#' (c2 = fmx_constraint(d2))
#' stopifnot(identical(c2, fmx_constraint_user(distname = 'GH', K = 2L, user = 'g2')))
#' fmx_constraint_brief(d2)
#' 
#' fmx_constraint_brief(fmx('norm', mean = c(0, 1)))
#' 
#' @name fmx_constraint
#' @export
fmx_constraint <- function(dist, distname = dist@distname, K = dim(dist@parM)[1L], parM = dist@parM) {
  switch(distname, norm = {
    return(integer())
  }, GH = {
    colID <- c('g', 'h')
    parM0 <- which((parM[, colID, drop = FALSE] == 0), arr.ind = TRUE)
    if (!length(parM0)) return(integer())
    gid <- which(parM[,'g'] == 0)
    gid1 <- 2L*K + gid # `g` parameters located at (2K+1L):(3K)
    hid <- which(parM[,'h'] == 0)
    hid1 <- 3L*K + hid # `h` parameters located at (3K+1L):(4K)
    ret <- c(gid1, hid1)
    attr(ret, 'user') <- paste0(colID[parM0[,'col']], parM0[,'row'])
    attr(ret, 'latex') <- paste0(colID[parM0[,'col']], '_{', parM0[,'row'], '}')
    attr(ret, 'gid') <- gid 
    attr(ret, 'hid') <- hid
    return(ret)
  }, stop('distname not supported yet:', sQuote(distname)))
}

#' @rdname fmx_constraint
#' @export
fmx_constraint_user <- function(distname, K, user) {
  switch(distname, norm = {
    return(integer())
  }, GH = {
    colID <- c('g', 'h')
    gid <- as.integer(gsub('^g', replacement = '', x = grep('^g', x = user, value = TRUE))) # len-0 compatible
    if (any(gid > K)) stop('having only ', K, ' components')
    gid1 <- 2L*K + gid # `g` parameters located at (2K+1L):(3K)
    hid <- as.integer(gsub('^h', replacement = '', x = grep('^h', x = user, value = TRUE))) # len-0 compatible
    if (any(hid > K)) stop('having only ', K, ' components')
    hid1 <- 3L*K + hid # `h` parameters located at (3K+1L):(4K)
    if (!length(ret <- c(gid1, hid1))) return(integer())
    attr(ret, 'user') <- user
    attr(ret, 'latex') <- c(if (length(gid)) paste0('g_{', gid, '}'), if (length(hid)) paste0('h_{', hid, '}'))
    attr(ret, 'gid') <- gid 
    attr(ret, 'hid') <- hid
    return(ret)
  }, stop('distname not supported yet:', sQuote(distname)))
}

#' @rdname fmx_constraint
#' @export
fmx_constraint_brief <- function(dist) {
  distname <- dist@distname
  K <- dim(dist@parM)[1L]
  distK <- paste0(distname, K)
  switch(distname, norm = return(distK), GH = {
    constr <- fmx_constraint(dist)
    usr <- attr(constr, which = 'user', exact = TRUE)
    if (!length(usr)) return(paste0('Unconstrained ', distK))
    if (identical(usr, c(t.default(outer(c('g', 'h'), seq_len(K), FUN = paste0))))) return(paste0('norm', K))
    latex <- attr(constr, which = 'latex', exact = TRUE)
    return(paste0(distK, ': $', paste(c(latex, '0'), collapse = '='), '$'))
  }, stop('unsupported distribution ', sQuote(distname)))
}




