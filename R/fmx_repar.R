




#' @title Reparameterization of \linkS4class{fmx} Object
#' 
#' @description 
#' To convert the parameters of \linkS4class{fmx} object into unrestricted parameters.
#' 
#' @param x \linkS4class{fmx} object
#' 
#' @param distname \link[base]{character} scalar, default `x@@distname`
#' 
#' @param pars \link[base]{numeric} \link[base]{matrix}, default `x@@pars`
#' 
#' @param K \link[base]{integer} scalar, default value from `x`
#' 
#' @param w \link[base]{numeric} \link[base]{vector}, default `x@@w`
#' 
#' @param ... additional parameters, not currently used
#' 
#' @details 
#' 
#' For the first parameter
#' \itemize{
#' \item {\eqn{A_1 \rightarrow A_1}}
#' \item {\eqn{A_2 \rightarrow A_1 + \exp(\log(d_1))}}
#' \item {\eqn{A_k \rightarrow A_1 + \exp(\log(d_1)) + \cdots + \exp(\log(d_{k-1}))}}
#' }
#' 
#' For mixing proportions to multinomial logits.
#' 
#' For `'norm'`: `sd -> log(sd)`
#' for `'GH'`: `B -> log(B), h -> log(h)`
#' 
#' @returns 
#' Function [fmx2dbl()] returns a \link[base]{numeric} \link[base]{vector}
#' 
#' @seealso [dbl2fmx()]
#' @export
fmx2dbl <- function(x, distname = x@distname, pars = x@pars, K = dim(pars)[1L], w = x@w, ...) { 
  # no longer used in compute intensive algorithms
  w_val <- qmlogis_first(w) # \code{K == 1L} will return \code{numeric(0)}
  if (!all(is.finite(w_val))) stop('NA or Inf in proportion indicated degenerated mixture (one or more component has 0% mixture proportion)')
  w_nm <- if (K == 1L) character() else paste0('logit', 2:K)
  pars[, id] <- log(pars[, (id <- dist_logtrans(distname))])
  argnm <- switch(distname, norm = c('mean', 'sdlog'), GH = c('A', 'Blog', 'g', 'hlog'), stop('write more'))
  if (K > 1L) {
    # when some log(d) is negative and has too great absolute value, 
    # exp(log(d)) is numerically 0, and `pars` will not be strictly increasing
    # in the method below we wont be able to retrieve the real log(d), instead get -Inf.
    pars[2:K, 1L] <- log(pars[2:K, 1L] - pars[1:(K-1L), 1L])
    locnm <- c(paste0(argnm[1L], 1L), paste0('log\u0394', seq_len(K)[-1L]))
  } else locnm <- paste0(argnm[1L], seq_len(K))
  out <- c(pars, w_val)
  names(out) <- c(locnm, paste0(rep(argnm[-1L], each = K), seq_len(K)), w_nm)
  return(out)
}



#' @title Inverse of [fmx2dbl()], for internal use
#' 
#' @description ..
#' 
#' @param x \link[base]{numeric} \link[base]{vector}, unrestricted parameters
#' 
#' @param K \link[base]{integer} scalar
#' 
#' @param distname \link[base]{character} scalar
#' 
#' @param ... additional parameters, not currently used
#' 
#' @details 
#' Only used in function [QLMDe()] and unexported function `qfmx_gr()`, not compute intensive
#' 
#' @returns 
#' Function [dbl2fmx()] returns a \link[base]{list} with two elements `$pars` and `$w`
#' 
#' @export
dbl2fmx <- function(x, K, distname, ...) {
  nx <- length(x)
  n_dist <- nx - (K - 1L) # K == 1L or not
  w <- if (K == 1L) 1 else unname(pmlogis_first(x[(n_dist + 1L):nx]))
  pm <- array(x[seq_len(n_dist)], dim = c(K, n_dist/K)) # not compute intensive..
  pm[,id] <- exp(pm[, (id <- dist_logtrans(distname)), drop = FALSE])
  if (K > 1L) pm[,1L] <- cumsum(c(pm[1L,1L], exp(pm[2:K,1L])))
  list(pars = pm, w = w) # much faster than ?base::cbind
}

