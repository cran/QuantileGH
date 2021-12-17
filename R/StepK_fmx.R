


#' @title Forward-Backward Selection of the Number of Components \eqn{K}
#' 
#' @description
#' 
#' The function \code{\link{StepK_fmx}} compares \eqn{gh}-parsimonious models with different number of components \eqn{K}
#' and selects the optimal model using the Vuong's closeness test.
#' 
#' @param object \code{'fmx_QLMDe'} object returned from \code{\link{QLMDe}}
#' 
#' @param by criterion, currently supporting 
#' \code{'logLik'} (via my \code{\link{LikRatio}}, likelihood ratio test, default and recommended), 
#' \code{'AIC'} (via \code{\link[stats]{AIC}}) and 
#' \code{'BIC'} (via \code{\link[stats]{BIC}}).
#' 
#' @param Kmax 'integer' value, the maximum number of component to be considered
#' 
#' @param silent 'logical' value, whether messages should be suppressed (default \code{FALSE})
#' 
#' @param ... additional parameters
#' 
#' @details 
#' 
#' The algorithm starts with selection \eqn{gh}-parsimonious model with user-specified initial number of components \eqn{K_0}.
#' Then the number of components is increased by \eqn{1} and the corresponding 
#' \eqn{gh}-parsimonious model is compared to the \eqn{gh}-parsimonious model with \eqn{K_0} components using the Vuong's closeness test. 
#' If \eqn{gh}-parsimonious model with \eqn{K_0} components is preferred then the algorithm is stopped if \eqn{K_0=1} 
#' or switches to backward selection if \eqn{K_0 > 1}.
#' If \eqn{gh}-parsimonious model with \eqn{K_0 + 1} components is preferred then
#' the algorithm is stopped if \eqn{K_0+1=Kmax} (prespecified maximum number of components)
#' or the next iteration of the algorithm is performed if \eqn{K_0+1<Kmax}.
#' The backward selection is performed only if \eqn{gh}-parsimonious model with \eqn{K_0} components 
#' is preferred to \eqn{gh}-parsimonious model with \eqn{K_0+1} components. 
#' Then \eqn{gh}-parsimonious model with \eqn{K_0-1} components is compared to
#' \eqn{gh}-parsimonious model with \eqn{K_0} components. 
#' If \eqn{gh}-parsimonious model with \eqn{K_0} components is preferred then the algorithm is stopped
#' and \eqn{gh}-parsimonious model with \eqn{K_0} components is optimal.
#' If \eqn{gh}-parsimonious model with \eqn{K_0-1} components is preferred then
#' the algorithm is stopped if \eqn{K_0-1=1} 
#' or the next iteration of the algorithm is performed if \eqn{K_0-1>1}.
#' 
#' @return 
#' 
#' \code{\link{StepK_fmx}} returns an \code{'fmx_QLMDe'} object, with attributes
#' \itemize{
#' \item{\code{anova}} {ANOVA table}
#' \item{\code{objF}} {value of the objective function (either the log-likelihood, AIC or BIC)}
#' }
#' 
#' @export
StepK_fmx <- function(object, by = c('logLik', 'AIC', 'BIC'), Kmax = 3L, silent = FALSE, ...) {
  by <- match.arg(by)
  K <- K_orig <- dim(object@parM)[1L]
  if (!silent) cat('Finding parsimonious model at K =', K, '\n')
  modelK <- model_orig <- Step_fmx(object, by = by, silent = silent) # comparison is to parsimonious model at original `K`
  objF <- objF_orig <- attr(model_orig, which = 'objF', exact = TRUE)[1L] # 'list'; objective function for parsimonious model at original `K`
  aod <- aod_orig <- attr(model_orig, which = 'anova', exact = TRUE)[1L, ] # # aod line for parsimonious model at original `K`
  
  compareK <- function(model0, model1) {
    # K0 > K1
    K0 <- dim(model0@parM)[1L]
    K1 <- dim(model1@parM)[1L]
    # return `TRUE` indicates selecting `K1` (K_less)
    if (K0 <= K1) stop('`model0` should have higher `K` than `model1`')
    #cat('larger model K0 = ', K0, '\n')
    #cat('smaller model K1 = ', K1, '\n')
    if ((distname <- model0@distname) != model1@distname) stop('`model0` has different `distname` than `model1`')
    objF0 <- attr(model0, which = 'objF', exact = TRUE)[1L]
    objF1 <- attr(model1, which = 'objF', exact = TRUE)[1L]
    switch(by, AIC = , BIC = {
      # smaller the better
      return(unlist(objF1) <= (unlist(objF0) + 1e-07)) # 'TRUE' for selecting `K1` (K_less)
    }, logLik = {
      lr_K <- LikRatio(models = c(objF0, objF1), # `model0` is more complicated than `model1`
                       type = switch(distname, norm = 'plain', GH = 'vuong', stop('unsupported ', sQuote(object@distname))), 
                       compare = 'first')
      if (inherits(lr_K, what = 'vuong')) {
        # see ?LikRatio: BIC-correction more prone to `model1` (simpler model)
        # 1-comp GH: up to 4; 2-comp normal: 2*2+1=5
        # 2-comp GH: up to 4*2+1=9; 3-comp normal: 3*2+2=8
        return(lr_K[2L, c('Decision_BIC')] != names(objF0))
        # if the more-prone-to-simpler-model-decision suggests `K0` (K_more), then choose K_more
        # otherwise, either '[tie]' (in df, i.e., number of parameters) or `K1` (K_less), then choose K_less
      } else if (inherits(lr_K, what = 'anova')) {
        return(lr_K[2L, length(lr_K)] > .05) # p-value of K_less; 'TRUE' for selecting `K1` (K_less)
      } else stop('should not come here')
      # `lr_K` is not returned, for now
    }, stop('unsupported ', sQuote(by)))
    return(FALSE)
  }
  
  ###################################
  # first, increase `K` (always compare to K-1, not `K_orig`)
  
  while ((K + 1L) <= Kmax) {
    tmp_old <- modelK
    obj_K <- QLMDe(object@data, distname = object@distname, data.name = object@data.name, K = K + 1L, p = object@p)
    if (!silent) cat('Finding parsimonious model at K =', K + 1L, '\n')
    tmp_new <- Step_fmx(obj_K, by = by, silent = silent, ...)
    #if (FALSE) {
    #  tmp_old <<- tmp_old
    #  tmp_new <<- tmp_new
    #  stop('here :)')
    #}
    if (compareK(model0 = tmp_new, model1 = tmp_old)) break # smaller model (i.e. `tmp_old`) is selected 
    K <- K + 1L
    if (!silent) cat('Increased K = ', K, ' is selected.\n', sep = '')
    modelK <- tmp_new
    aod <- rbind.data.frame(aod, attr(modelK, which = 'anova', exact = TRUE)[1L, ]) # still 'anova'
    objF <- c(objF, attr(modelK, which = 'objF', exact = TRUE)[1L]) # 'list'
  }
  
  ###################################
  # then, decrease `K`
  
  if (dim(modelK@parM)[1L] == K_orig) { # only do `descreasing` if no `increasing`
    # always compare to `K_orig`, not K-1
    while (K > 1L) {
      obj_K <- QLMDe(object@data, distname = object@distname, data.name = object@data.name, K = K - 1L, p = object@p)
      if (!silent) cat('Finding parsimonious model at K =', K - 1L, '\n')
      tmpK <- Step_fmx(obj_K, by = by, silent = silent, ...)
      if (!compareK(model0 = model_orig, model1 = tmpK)) break # larger model (i.e. `model_orig` is retained)
      K <- K - 1L
      if (!silent) cat('Reduced K = ', K, ' is selected.\n', sep = '')
      modelK <- tmpK
      aod <- rbind.data.frame(aod, attr(modelK, which = 'anova', exact = TRUE)[1L, ]) # still 'anova'
      objF <- c(objF, attr(modelK, which = 'objF', exact = TRUE)[1L]) # 'list'
    }
  }
  
  attr(aod, 'heading') <- 'Stepwise #-of-Comp Selection (Parsimonious per #-of-Comp)'
  attr(modelK, 'anova') <- aod
  attr(modelK, 'objF') <- objF
  if (dim(modelK@parM)[1L] != K_orig) attr(modelK, 'orig_K') <- model_orig # parsimonious model at original K
  return(modelK)
  
}






if (FALSE) {
  (d = fmx('norm', mean = c(1, 4, 8), w = c(3, 3, 4)))
  x = rfmx(n = 1e3L, dist = d)
  gghist(x)
  y1 = QLMDe(x, distname = 'norm', K = 1L)
  StepK_fmx(y1, Kmax = 3L)
  
  (d = fmx('GH', A = c(1, 4, 8), g = c(.2, 0, -.3), h = c(.2, .2, .2), w = c(3, 3, 4)))
  x = rfmx(n = 1e3L, dist = d)
  gghist(x)
  (y1 = QLMDe(x, distname = 'GH', K = 1L))
  StepK_fmx(y1)
}





