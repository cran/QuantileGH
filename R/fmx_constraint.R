

#' @title Parameter Constraint(s) of Mixture Distribution
#' 
#' @description 
#' 
#' Determine the parameter constraint(s) of a finite mixture distribution \linkS4class{fmx}, 
#' either by the value of parameters of such mixture distribution, 
#' or by a user-specified string. 
#' 
#' @param dist (optional) \linkS4class{fmx} object
#' 
#' @param distname \link[base]{character} scalar, name of distribution (see \linkS4class{fmx}),
#' default value determined by `dist`
#' 
#' @param K \link[base]{integer} scalar, number of components,
#' default value determined by `dist`
#' 
#' @param pars \link[base]{double} \link[base]{matrix}, 
#' distribution parameters of a finite mixture distribution (see \linkS4class{fmx}),
#' default value determined by `dist`
#' 
#' @returns 
#' 
#' \link{fmx_constraint} returns the indices of internal parameters 
#' (only applicable to Tukey's \eqn{g}-&-\eqn{h} mixture distribution, yet) to be constrained, 
#' based on the input \linkS4class{fmx} object `dist`.
#' 
#' @examples 
#' (d0 = fmx('GH', A = c(1,4), g = c(.2,.1), h = c(.05,.1), w = c(1,1)))
#' (c0 = fmx_constraint(d0))
#' stopifnot(identical(c0, user_constraint(character(), distname = 'GH', K = 2L)))
#' 
#' (d1 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(0,.1), w = c(1,1)))
#' (c1 = fmx_constraint(d1))
#' stopifnot(identical(c1, user_constraint(c('g2', 'h1'), distname = 'GH', K = 2L)))
#' 
#' (d2 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(.15,.1), w = c(1,1)))
#' (c2 = fmx_constraint(d2))
#' stopifnot(identical(c2, user_constraint('g2', distname = 'GH', K = 2L)))
#' 
#' @name fmx_constraint
#' @export
fmx_constraint <- function(dist, distname = dist@distname, K = dim(dist@pars)[1L], pars = dist@pars) {
  switch(distname, GH = {
    colID <- c('g', 'h')
    pars0 <- which((pars[, 3:4, drop = FALSE] == 0), arr.ind = TRUE)
    if (!length(pars0)) return(integer())
    gid <- which(pars[,3L] == 0)
    gid1 <- 2L*K + gid # `g` parameters located at (2K+1L):(3K)
    hid <- which(pars[,4L] == 0)
    hid1 <- 3L*K + hid # `h` parameters located at (3K+1L):(4K)
    ret <- c(gid1, hid1)
    attr(ret, which = 'user') <- paste0(colID[pars0[,'col']], pars0[,'row'])
    attr(ret, which = 'tex') <- paste0(colID[pars0[,'col']], '_{', pars0[,'row'], '}')
    attr(ret, which = 'gid') <- gid 
    attr(ret, which = 'hid') <- hid
    return(ret)
  }, return(integer()))
}





#' @title Formalize User-Specified Constraint of \linkS4class{fmx} Object
#' 
#' @description 
#' Formalize user-specified constraint of \linkS4class{fmx} object
#' 
#' @param x \link[base]{character} \link[base]{vector}, constraint(s) to be imposed.  
#' For example, for a two-component Tukey's \eqn{g}-&-\eqn{h}
#' mixture, `c('g1', 'h2')` indicates \eqn{g_1=h_2=0} given \eqn{A_1 < A_2}, i.e., the 
#' \eqn{g}-parameter for the first component (with smaller location value)
#' and the \eqn{h}-parameter for the second component (with larger mean value) are to be constrained as 0.
#' 
#' @param distname \link[base]{character} scalar, name of distribution
#' 
#' @param K \link[base]{integer} scalar, number of components
#' 
#' @returns 
#' 
#' \link{user_constraint} returns the indices of internal parameters 
#' (only applicable to Tukey's \eqn{g}-&-\eqn{h} mixture distribution, yet) to be constrained, 
#' based on the type of distribution `distname`, number of components `K`
#' and a user-specified string (e.g., `c('g2', 'h1')`).
#' 
#' @examples 
#' (d0 = fmx('GH', A = c(1,4), g = c(.2,.1), h = c(.05,.1), w = c(1,1)))
#' (c0 = fmx_constraint(d0))
#' stopifnot(identical(c0, user_constraint(distname = 'GH', K = 2L, x = character())))
#' 
#' (d1 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(0,.1), w = c(1,1)))
#' (c1 = fmx_constraint(d1))
#' stopifnot(identical(c1, user_constraint(distname = 'GH', K = 2L, x = c('g2', 'h1'))))
#' 
#' (d2 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(.15,.1), w = c(1,1)))
#' (c2 = fmx_constraint(d2))
#' stopifnot(identical(c2, user_constraint(distname = 'GH', K = 2L, x = 'g2')))
#' 
#' @export
user_constraint <- function(x, distname, K) {
  switch(distname, GH = {
    colID <- c('g', 'h')
    gid <- as.integer(gsub('^g', replacement = '', x = grep('^g', x = x, value = TRUE))) # len-0 compatible
    if (any(gid > K)) stop('having only ', K, ' components')
    gid1 <- 2L*K + gid # `g` parameters located at (2K+1L):(3K)
    hid <- as.integer(gsub('^h', replacement = '', x = grep('^h', x = x, value = TRUE))) # len-0 compatible
    if (any(hid > K)) stop('having only ', K, ' components')
    hid1 <- 3L*K + hid # `h` parameters located at (3K+1L):(4K)
    if (!length(ret <- c(gid1, hid1))) return(integer())
    attr(ret, which = 'user') <- x
    attr(ret, which = 'tex') <- c(if (length(gid)) paste0('g_{', gid, '}'), if (length(hid)) paste0('h_{', hid, '}'))
    attr(ret, which = 'gid') <- gid 
    attr(ret, which = 'hid') <- hid
    return(ret)
  }, return(integer()))
}





#' @title TeX Label (of Parameter Constraint(s)) of \linkS4class{fmx} Object
#' 
#' @description 
#' 
#' Create TeX label of (parameter constraint(s)) of \linkS4class{fmx} object
#' 
#' @param dist \linkS4class{fmx} object
#' 
#' @param print_K \link[base]{logical} scalar, whether to print the number of components \eqn{K}.
#' Default `FALSE`.
#' 
#' @returns 
#' 
#' Function [getTeX()] returns a \link[base]{character} scalar 
#' (of TeX expression) of the constraint, 
#' primarily intended for end-users in plots.
#' 
#' 
#' @examples 
#' (d0 = fmx('GH', A = c(1,4), g = c(.2,.1), h = c(.05,.1), w = c(1,1)))
#' getTeX(d0)
#' 
#' (d1 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(0,.1), w = c(1,1)))
#' getTeX(d1)
#' 
#' (d2 = fmx('GH', A = c(1,4), g = c(.2,0), h = c(.15,.1), w = c(1,1)))
#' getTeX(d2)
#' 
#' @export
getTeX <- function(dist, print_K = FALSE) {
  if (!inherits(dist, what = 'fmx')) stop('do not allow')
  distname <- dist@distname
  K <- dim(dist@pars)[1L]
  distK <- paste0(distname, K)
  switch(distname, GH = {
    constr <- fmx_constraint(dist)
    usr <- attr(constr, which = 'user', exact = TRUE)
    if (!length(usr)) return(paste0('Full ', distK))
    if (identical(usr, c(t.default(outer(c('g', 'h'), seq_len(K), FUN = paste0))))) return(paste0('norm', K))
    latex <- attr(constr, which = 'tex', exact = TRUE)
    ret <- paste(c(latex, '0'), collapse = '=')
    if (print_K) return(paste0(distK, ': $', ret, '$')) #return(paste0('$', ret, ', K=', K, '$'))
    return(paste0('$', ret, '$'))
  }, { # all else
    nm <- switch(distname, 
                 norm = 'Normal',
                 sn = 'Skew Normal',
                 st = 'Skew $t$', 
                 genpois1 = 'Generalized Poisson', # ?SIBERG::fitGP
                 nbinom = 'Negative Binomial', # SIBERG::fitNB
                 stop())
    if (print_K) sprintf('%d-comp. %s', K, nm) else nm
  })
}




#' @title Number of Parameters of \linkS4class{fmx} Object
#' 
#' @description ..
#' 
#' @param dist \linkS4class{fmx} object
#' 
#' @details 
#' Also the degree-of-freedom in \link[stats]{logLik},
#' as well as \code{stats:::AIC.logLik} and \code{stats:::BIC.logLik}
#' 
#' @returns 
#' \link{npar.fmx} returns an \link[base]{integer} scalar.
#' 
#' @export
npar.fmx <- function(dist) {
  # https://en.wikipedia.org/wiki/Akaike_information_criterion
  # ?stats::logLik
  # ?stats:::AIC.default
  # attr(, 'df') is the number of (estimated) parameters in the model.
  # I write this function so that I do not have to do \link{dfmx} if not needed
  dm <- dim(dist@pars)
  ret <- (dm[2L] + 1L) * dm[1L] - 1L - length(attr(fmx_constraint(dist), which = 'user', exact = TRUE))
  if (dist@distname == 'st') {
    nu <- dist@pars[,'nu']
    if (length(unique.default(nu)) != 1L) stop('\\pkg{mixsmsn} update to enable multiple `nu`?')
    ret <- ret - (length(nu) - 1L)
  }
  return(ret)
}




