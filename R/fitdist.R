
#' @title Log-Likelihood of \link[fitdistrplus]{fitdist} Object
#' 
#' @description ..
#' 
#' @param object \link[fitdistrplus]{fitdist} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' Output of \link[fitdistrplus]{fitdist} has elements `$loglik`, `$aic` and `$bic`, 
#' but they are simply \link[base]{numeric} scalars.
#' `fitdistrplus:::logLik.fitdist` simply returns these elements.
#' 
#' Function [logLik.fitdist()] returns a \link[stats]{logLik} object, which 
#' could be further used by \link[stats]{AIC} and \link[stats]{BIC}.
#' 
#' (I have written to the authors)
#' 
#' @returns 
#' Function [logLik.fitdist()] returns a \link[stats]{logLik} object
#' 
#' @importFrom stats logLik
#' @export logLik.fitdist
#' @export
logLik.fitdist <- function(object, ...) {
  ret <- object[['loglik']]
  attr(ret, which = 'nobs') <- object[['n']]
  attr(ret, which = 'df') <- length(object[['estimate']])
  class(ret) <- 'logLik'
  return(ret)
}


#' @title Number of Observations in \link[fitdistrplus]{fitdist} Object
#' 
#' @description ..
#' 
#' @param object \link[fitdistrplus]{fitdist} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @returns 
#' Function [nobs.fitdist()] returns an \link[base]{integer} scalar
#' 
#' @importFrom stats nobs
#' @export nobs.fitdist
#' @export
nobs.fitdist <- function(object, ...) object[['n']]



# ?fitdistrplus:::coef.fitdist # slow
