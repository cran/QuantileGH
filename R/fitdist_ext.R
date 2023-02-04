
#' @title Log-Likelihood of \link[fitdistrplus]{fitdist} Object
#' 
#' @description ..
#' 
#' @param object \link[fitdistrplus]{fitdist} object
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' Output of \link[fitdistrplus]{fitdist} has element \code{$loglik} 
#' (as well as \code{$aic} and \code{$bic}), 
#' but it is simply a \link[base]{numeric} scalar.
#' \code{fitdistrplus:::logLik.fitdist} simply returns this element.
#' 
#' \link{logLik.fitdist} returns a \link[stats]{logLik} object, which 
#' could be further used by \link[stats]{AIC} and \link[stats]{BIC}.
#' 
#' (I have written to the authors)
#' 
#' @return 
#' \link{logLik.fitdist} returns a \link[stats]{logLik} object
#' 
#' @seealso \link[fitdistrplus]{fitdist}
#' 
#' @importFrom stats logLik
#' 
#' @export
logLik.fitdist <- function(object, ...) {
  n <- length(data <- object[['data']])
  if (!n) stop('Re-run ?fitdistrplus::fitdist with `keepdata = TRUE`')
  ret <- object[['loglik']]
  attr(ret, which = 'nobs') <- n
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
#' @return 
#' \link{nobs.fitdist} returns an \link[base]{integer} scalar
#' 
#' @seealso \link[fitdistrplus]{fitdist}
#' 
#' @importFrom stats nobs
#' 
#' @export
nobs.fitdist <- function(object, ...) object[['n']]



# ?fitdistrplus:::coef.fitdist # slow
