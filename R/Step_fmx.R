

# S3 method for ?stats::add1
#' @export
add1.fmx_QLMDe <- function(object, scope, ...) {
  stop('not written yet')
}


# S3 method for ?stats::drop1
# do \strong{not} write as S3 of ?MASS::dropterm; there's no \strong{term} for 'fmx_QLMDe' object
#' @export
drop1.fmx_QLMDe <- function(object, scope, silent = FALSE, ...) { # `...` not used
  K <- dim(object@parM)[1L]
  nms <- paste0(K, '-comp ', object@distname, ': ', vapply(scope, FUN = \(i) paste(c(i,'0'), collapse = '='), FUN.VALUE = ''))
  # input `object` is \strong{not} included in the elements of return
  lapply(setNames(seq_along(scope), nm = nms), FUN = \(i) {
    if (!silent) cat(nms[i], '.. ')
    ret <- QLMDe(object@data, data.name = object@data.name, distname = object@distname, K = K, p = object@p, constraint = scope[[i]])
    if (!silent) cat('done! (stats::optim iter. ', ret@optim$counts['function'], ')\n', sep = '')
    return(ret)
  })
}



#' @title Backward Selection \eqn{gh}-parsimonious Model with Fixed Number of Components
#' 
#' @description 
#' 
#' The function \code{\link{Step_fmx}} selects a \eqn{gh}-parsimonious model with \eqn{g} and/or \eqn{h} parameters equal to zero 
#' for all or some of the mixture components conditionally on fixed number of components \eqn{K}.
#' 
#' @param object \code{'fmx_QLMDe'} object returned from \code{\link{QLMDe}}
#' 
#' @param by criterion, currently supporting 
#' \code{'logLik'} (via my \code{\link{LikRatio}}, likelihood ratio test, default and recommended), 
#' \code{'AIC'} (via \code{\link[stats]{AIC}}) and 
#' \code{'BIC'} (via \code{\link[stats]{BIC}}).
#' 
#' @param silent 'logical' value, whether messages should be suppressed (default \code{FALSE})
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
#' \code{\link{Step_fmx}} returns an \code{'fmx_QLMDe'} object, with attributes
#' \itemize{
#' \item{\code{anova}} {ANOVA table}
#' \item{\code{objF}} {value of the objective function (either the log-likelihood, AIC or BIC)}
#' }
#' 
#' @export
Step_fmx <- function(object, by = c('logLik', 'AIC', 'BIC'), silent = FALSE, ...) {
  
  # 'fmx_fit' is for back-compatibility for simulation results
  # in future may add 'fmx_newMethod'
  if (!inherits(object, what = c('fmx_QLMDe', 'fmx_fit'))) stop('input must be \'fmx_QLMDe\'')
  by <- match.arg(by)
  
  K <- dim(object@parM)[1L]
  candpar <- switch(object@distname, GH = { # candidate parameter(s) to be eliminated
    c(outer(c('g', 'h'), 1:K, FUN = paste0))
  }, norm = character(), stop('unsupported ', sQuote(object@distname)))
  
  fit <- object
  fit_constr <- attr(fmx_constraint(object), which = 'user', exact = TRUE)
  
  if (length(fit_constr)) { # in future
    # then add1
  }
  
  npar <- npar_fmx(object)
  model_name <- paste0(K, '-comp ', switch(object@distname, GH = {
    paste0(object@distname, ': ', if (!length(fit_constr)) 'Unconstraint' else paste(c(fit_constr,0), collapse = '='))
  }, norm = 'Normal'))
  fit_orig <- setNames(list(object), nm = model_name)
  objF <- objF_orig <- setNames(list(eval(call(by, object))), nm = model_name)
  
  # comparison is to original `object`
  while (length(freepar <- setdiff(candpar, fit_constr))) { # still have `freepar` to be eliminated
    
    scope <- lapply(freepar, FUN = \(i) candpar[candpar %in% c(i, fit_constr)]) # preserve conventional order of parameter(s)
    tmp <- c(fit_orig, drop1(object = fit, scope = scope, silent = silent)) # ?drop1.fmx_QLMDe
    n_objF <- c(objF_orig, lapply(tmp[-1L], FUN = \(i) {
      tryCatch(eval(call(by, i)), error = as.null.default) # mal-fit, very very rare
    })) # new objF
    if (any(id <- (lengths(n_objF, use.names = FALSE) == 0L))) {
      cat('.. mal-fit:', paste(sQuote(names(n_objF)[id]), collapse = ', '), '(at %\'s of original fit)\n')
      n_objF <- n_objF[!id]
    }
    
    switch(by, logLik = {
      lr <- LikRatio(models = n_objF, type = 'plain', compare = 'first') # print(lr) # ?stats:::print.anova
      pval <- lr[[length(lr)]] # p-value must on last column (see my ?.pval.anova); 1st element of p-value always NA_real_
      o1 <- which.max(pval) # 1st element NA_real_ omitted automatically
      if (pval[o1] < .05) break # original `object` is significantly different from reduced model(s)
    }, AIC =, BIC = {
      o1 <- order(unlist(n_objF, use.names = FALSE))[1L]
      if (o1 == 1L) break # original `object` has smallest AIC/BIC
    }, stop('unsupported ', sQuote(by)))
    
    fit <- tmp[[o1]]
    model_name <- c(model_name, names(tmp)[o1])
    npar <- c(npar, npar_fmx(fit))
    objF <- c(objF, n_objF[o1]) # 'list'
    fit_constr <- attr(fmx_constraint(fit), which = 'user', exact = TRUE)
    
  }
  
  aod <- data.frame(
    '# Parameter' = rev.default(npar),
    objF = rev.default(unlist(objF, use.names = FALSE)), 
    check.names = FALSE, row.names = rev.default(model_name))
  names(aod)[2L] <- by 
  attr(aod, 'heading') <- 'Stepwise Parameter Elimination (Fixed # of Comp.)'
  class(aod) <- c('anova', 'data.frame')
  if (!silent) print(aod)
  attr(fit, 'anova') <- aod
  attr(fit, 'objF') <- rev.default(objF)
  return(fit)
}




