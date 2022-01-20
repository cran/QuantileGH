


#' @title Add or Drop All Possible Parameters of \code{\linkS4class{fmx_QLMDe}} Object
#' 
#' @description 
#' 
#' Compute all the single terms in the \code{scope} argument that can be added to or dropped from the model, 
#' fit those models and compute a table of the changes in fit.
#' 
#' These are additional S3 methods of \code{\link[stats]{add1}} and \code{\link[stats]{drop1}}.
#' 
#' @param object \code{\linkS4class{fmx_QLMDe}} object
#' 
#' @param scope a \code{\link[base]{list}} of \code{\link[base]{character}} vectors to denote one or more constraints
#' 
#' @param test \code{\link[base]{character}}, either \code{'logLik'} (default), \code{'AIC'} or \code{'BIC'}
#' 
#' @param trace \code{\link[base]{logical}} scalar, whether to print out progress reports (default \code{FALSE})
#' 
#' @param ... place holder to match the S3 generic \code{\link[stats]{drop1}}, currently not in use
#' 
#' @details 
#' 
#' Do \strong{not} write as S3 method of \code{\link[MASS]{dropterm}}; there's no \strong{term} for \code{\linkS4class{fmx_QLMDe}} object.
#' 
#' @return
#' 
#' \code{\link{drop1.fmx_QLMDe}} returns an \code{\link[stats]{anova}} table with additional attributes
#' \itemize{
#' \item{\code{models}} {a \code{\link[base]{list}} of \code{\linkS4class{fmx_QLMDe}} objects}
#' \item{\code{objF}} {a \code{\link[base]{list}} of objective functions (depends on \code{test})}
#' \item{\code{o1}} {the location of the optimal models by \code{test}.  If the original model is optimal, this value is \code{integer()}}
#' }
#'
#' \code{\link{add1.fmx_QLMDe}} will be added in the next release.
#' 
#' @seealso \code{\link[stats]{add1}}, \code{\link[stats]{drop1}}.
#' 
#' @name drop1_fmx
#' @export
drop1.fmx_QLMDe <- function(object, scope, test = c('logLik', 'AIC', 'BIC'), trace = TRUE, ...) { # `...` not used
  
  K <- dim(object@parM)[1L]
  
  y0 <- lapply(scope, FUN = function(iscope) { # (iscope = scope[[1L]])
    if (trace) cat(paste0(object@distname, K), paste(c(iscope, '0'), collapse = '='), '.. ')
    ret <- QLMDe(object@data, data.name = object@data.name, distname = object@distname, K = K, p = object@p, constraint = iscope)
    if (trace) cat('done!\n', sep = '')
    return(ret)
  })
  y <- c(list(object), y0) # input `object` \strong{is} included in the elements of return
  names(y) <- vapply(y, FUN = fmx_constraint_brief, FUN.VALUE = '', USE.NAMES = FALSE)
  
  lr <- LikRatio(models = y, type = 'plain', compare = 'first') # print(lr) # ?stats:::print.anova
  attr(lr, which = 'models') <- y
  
  test <- match.arg(test)
  attr(lr, which = 'o1') <- switch(test, logLik = {
    pval <- lr[[length(lr)]] # p-value must on last column (see my ?.pval.anova); 1st element of p-value always NA_real_
    o1 <- which.max(pval) # 1st element NA_real_ omitted automatically
    if (pval[o1] > .05) o1 else integer() # original `object` is significantly different from reduced model(s)
  }, AIC =, BIC = {
    o1 <- order(lr[[test]])[1L]
    if (o1 != 1L) o1 else integer() # original `object` has smallest AIC/BIC
  })
  
  return(lr)
  
}

#' @rdname drop1_fmx
#' @export
add1.fmx_QLMDe <- function(object, scope, ...) {
  stop('not written yet')
}




#' @title Backward Selection \eqn{gh}-parsimonious Model with Fixed Number of Components
#' 
#' @description 
#' 
#' \code{\link{Step_fmx}} selects a \eqn{gh}-parsimonious model with \eqn{g} and/or \eqn{h} parameters equal to zero 
#' for all or some of the mixture components conditionally on fixed number of components \eqn{K}.
#' 
#' @param object \code{\linkS4class{fmx_QLMDe}} object returned from \code{\link{QLMDe}}
#' 
#' @param test \code{\link[base]{character}} value indicating the criterion to be used, with options
#' \code{'logLik'} (via likelihood ratio test \code{\link{LikRatio}}, default and recommended), 
#' \code{'AIC'} (via Akaike's information criterion \code{\link[stats]{AIC}}) and 
#' \code{'BIC'} (via Bayesian information criterion \code{\link[stats]{BIC}}).
#' 
#' @param ... additional parameters
#' 
#' @details 
#' 
#' The algorithm starts with quantile least Mahalanobis distance estimates (\code{\link{QLMDe}}) 
#' of either the full mixture of Tukey \eqn{g}-&-\eqn{h} distributions model, or
#' a constrained model (i.e., some \eqn{g} and/or \eqn{h} parameters equal to zero according to the user input).
#' Next, each of the non-zero \eqn{g} and/or \eqn{h} parameters is tested using the likelihood ratio test.
#' If all tested \eqn{g} and/or \eqn{h} parameters are significantly different from zero at the level 0.05
#' the algorithm is stopped and the initial model is considered \eqn{gh}-parsimonious.
#' Otherwise, the \eqn{g} or \eqn{h} parameter with the largest p-value is constrained to zero 
#' for the next iteration of the algorithm.
#' 
#' The algorithm iterates until only significantly-different-from-zero \eqn{g} and \eqn{h} parameters 
#' are retained, which corresponds to \eqn{gh}-parsimonious Tukey's \eqn{g}-&-\eqn{h} mixture model.
#' 
#' @return 
#' 
#' \code{\link{Step_fmx}} returns an \code{\linkS4class{fmx_QLMDe}} object, with attributes
#' \itemize{
#' \item{\code{anova}} {ANOVA table}
#' \item{\code{objF}} {value of the objective function (either the log-likelihood, AIC or BIC)}
#' }
#' 
#' @export
Step_fmx <- function(object, test = c('logLik', 'AIC', 'BIC'), ...) {
  
  # in future may add 'fmx_newMethod'
  if (!inherits(object, what = c('fmx_QLMDe'))) stop('input must be \'fmx_QLMDe\'')
  test <- match.arg(test)
  
  K <- dim(object@parM)[1L]
  candpar <- switch(object@distname, GH = { # candidate parameter(s) to be eliminated
    c(outer(c('g', 'h'), 1:K, FUN = paste0))
  }, norm = character(), stop('unsupported ', sQuote(object@distname)))
  
  fit <- object
  fit_constr <- attr(fmx_constraint(fit), which = 'user', exact = TRUE)
  
  if (length(fit_constr)) { # in future
    # then add1
  }
  
  npar <- npar_fmx(fit)
  model_name <- fmx_constraint_brief(fit)
  #fit_orig <- setNames(list(fit), nm = model_name)
  objF <- setNames(list(eval(call(test, fit))), nm = model_name)
  
  # comparison vs original `object`
  while (length(freepar <- setdiff(candpar, fit_constr))) { # still have `freepar` to be eliminated
    scope <- lapply(freepar, FUN = function(i) candpar[candpar %in% c(i, fit_constr)]) # preserve conventional order of parameter(s)
    
    lr <- drop1.fmx_QLMDe(object = object, scope = scope, test = test, ...)
    if (!length(o1 <- attr(lr, which = 'o1', exact = TRUE))) break # original `object` is best by `test`
    
    imods <- attr(lr, which = 'models', exact = TRUE) 
    iobjF <- attr(lr, which = test, exact = TRUE)
    
    fit <- imods[[o1]]
    model_name <- c(model_name, names(imods)[o1])
    npar <- c(npar, npar_fmx(fit))
    objF <- c(objF, iobjF[o1]) # 'list'
    fit_constr <- attr(fmx_constraint(fit), which = 'user', exact = TRUE)
    
  }
  
  aod <- data.frame(
    '# Parameter' = rev.default(npar),
    objF = rev.default(unlist(objF, use.names = FALSE)), 
    check.names = FALSE, row.names = rev.default(model_name))
  names(aod)[2L] <- test 
  attr(aod, 'heading') <- 'Stepwise Parameter Elimination (Fixed # of Comp.)'
  class(aod) <- c('anova', 'data.frame')
  attr(fit, 'anova') <- aod
  attr(fit, 'objF') <- rev.default(objF)
  return(fit)
}




