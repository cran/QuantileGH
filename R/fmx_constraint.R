

#' @title Parameter Constraint(s) of Mixture Distribution
#' 
#' @description 
#' 
#' Determine the parameter constraint(s) of a finite mixture distribution, either 
#' by the value of parameters of such mixture distribution, or by a user-specified string. 
#' 
#' @param dist an \linkS4class{fmx} object, can be missing
#' 
#' @param distname \link[base]{character} scalar, name of distribution
#' 
#' @param K \link[base]{integer} scalar, number of components
#' 
#' @param parM \link[base]{double} \link[base]{matrix}, 
#' distribution parameters of a finite mixture distribution as slot \code{@@parM} of \linkS4class{fmx} object.
#' 
#' @param user \link[base]{character} vector, constraint(s) to be 
#' imposed for an \linkS4class{fmx} object.  For example, for a two-component Tukey's \eqn{g}-&-\eqn{h}
#' mixture, \code{user = c('g2', 'h1')} indicates the \eqn{g}-parameter for the first component (with smaller mean value)
#' and the \eqn{h}-parameter for the second component (with larger mean value) are to be constrained, i.e., \eqn{g_2=h_1=0}.
#' 
#' @return 
#' 
#' \link{fmx_constraint} returns the indexes of internal parameters 
#' (only applicable to Tukey's \eqn{g}-&-\eqn{h} mixture distribution, yet) to be constrained, 
#' based on the input \linkS4class{fmx} object \code{dist}.
#' 
#' \link{fmx_constraint_user} returns the indexes of internal parameters 
#' (only applicable to Tukey's \eqn{g}-&-\eqn{h} mixture distribution, yet) to be constrained, 
#' based on the type of distribution (\code{distname}), number of components (\code{K}) 
#' and a user-specified string (e.g., \code{c('g2', 'h1')}).
#' 
#' \link{fmx_constraint_brief} returns a \link[base]{character} scalar (of LaTeX expression) of the constraint, 
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
    #parM0 <- which((parM[, colID, drop = FALSE] == 0), arr.ind = TRUE)
    parM0 <- which((parM[, 3:4, drop = FALSE] == 0), arr.ind = TRUE)
    if (!length(parM0)) return(integer())
    #gid <- which(parM[,'g'] == 0)
    gid <- which(parM[,3L] == 0)
    gid1 <- 2L*K + gid # `g` parameters located at (2K+1L):(3K)
    #hid <- which(parM[,'h'] == 0)
    hid <- which(parM[,4L] == 0)
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




