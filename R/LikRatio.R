

# '\u0394' is lower case delta # unicode creates a NOTE in R CMD check

#' @title Likelihood Ratio Test for General Models
#' 
#' @description 
#' 
#' Likelihood ratio test for models fitted by R.
#' 
#' @param dots a \link[base]{list} of regression models, or a \link[base]{list} of \link[stats]{logLik} objects.
#' 
#' @param type \link[base]{character} scalar, ordinary likelihood ratio test (\code{'plain'}, default)
#' or Vuong's closeness test for non-nested models (\code{'vuong'}).
#' 
#' @param compare type of comparison between the models, sequentially (\code{'seq'}, default) or 
#' all models versus the first model (\code{'first'})
#' 
#' @param ... additional arguments of \link[stats]{logLik} function(s)
#' 
#' @seealso \link[lmtest]{lrtest.default}
#' 
#' @return 
#' 
#' \link{LikRatio} returns an \link[stats:anova]{ANOVA} table for likelihood ratios test, 
#' or a \code{'vuong'} object for Vuong's test.
#' 
#' @examples 
#' # no examples for now
#' 
#' @references 
#' Vuong's closeness test, \doi{10.2307/1912557}.
#' 
#' @export
LikRatio <- function(dots, type = c('plain', 'vuong'), compare = c('seq', 'first'), ...) {
  
  if (!is.list(dots) || any(!lengths(dots, use.names = FALSE))) stop('must manually remove NULL elements in `dots`')
  if ((n <- length(dots)) < 2L) stop('Must provide 2 or more models in `dots`')
  nm <- names(dots)
  if (!length(nm) || anyNA(nm) || !all(nzchar(nm))) stop('`dots` must be fully named')
  
  .logL <- lapply(dots, FUN = logLik, ...) # ?stats:::logLik.logLik
  .AIC <- lapply(.logL, FUN = AIC) # ?stats:::AIC.logLik
  .BIC <- lapply(.logL, FUN = BIC) # ?stats:::BIC.logLik
  
  Ns <- vapply(.logL, FUN = attr, which = 'nobs', exact = TRUE, FUN.VALUE = 0L)
  if (!all(duplicated.default(Ns)[-1L])) stop('sample size must be same')
  N <- Ns[1L] # needed in 'vuong'
  
  id0 <- switch(match.arg(compare), seq = 1:(n - 1L), first = 1L)
  id1 <- 2:n
  
  df <- vapply(.logL, FUN = attr, which = 'df', exact = TRUE, FUN.VALUE = 0)
  logL0 <- unlist(.logL, use.names = FALSE)
  
  d_df <- df[id1] - df[id0] # delta-df
  d_logL <- logL0[id1] - logL0[id0] # delta-log-likelihood
  
  out0 <- data.frame( # colnames from ?lmtest::lrtest.default
    Df = df, 
    'delta-Df' = c(NA_real_, d_df),
    logLik = logL0, 
    AIC = unlist(.AIC, use.names = FALSE), 
    BIC = unlist(.BIC, use.names = FALSE), 
    row.names = c(nm[1L], paste0(nm[id1], ' (vs. ', nm[id0], ')')),
    check.names = FALSE)
  
  switch(match.arg(type), plain = {
    
    chisq <- 2 * abs(d_logL)
    pval <- ifelse((.d_df <- round(abs(d_df))) == 0, yes = NA_real_, no = pchisq(chisq, df = .d_df, lower.tail = FALSE))
    out <- data.frame(
      out0,
      '2*|delta-logLik|' = c(NA_real_, chisq),
      'Pr(>Chisq)' = c(NA_real_, pval),
      check.names = FALSE)
    attr(out, which = 'heading') <- 'Likelihood ratio test'
    class(out) <- c('anova', 'data.frame') # can invoke ?stats:::print.anova
    
  }, vuong = {
    
    # \url{https://en.wikipedia.org/wiki/Vuong%27s_closeness_test}
    # Vuong (1989) https://authors.library.caltech.edu/81424/1/sswp605.pdf
    # Note that the numerator is the difference between \strong{log}-likelihoods.
    logl <- lapply(.logL, FUN = attr, which = 'logl', exact = TRUE)
    if (!all(lengths(logl, use.names = FALSE) == N)) stop('`logLik.*` must return pointwise log-likelihood')
    
    V_denom <- sqrt(N) * mapply(FUN = function(e1, e2) {
      # var(e1 - e2)
      mad(e1 - e2)^2 # use stats::mad()^2 to replace stats::var for robustness; important!!
    } , e1 = logl[id1], e2 = logl[id0], SIMPLIFY = TRUE)
    
    # both AIC- and BIC- correction are mentioned in Vuong (1989) Page 27, bottom
    z_AIC <- (d_logL - d_df) / V_denom # AIC-corrected
    z_BIC <- (d_logL - d_df * log(N) / 2) / V_denom # BIC-corrected
    # presume that `model0` is more complicated than `model1` (i.e., in ?StepK_fmx), then d_df < 0
    # since log(N)/2 almost always > 1 (i.e., N > 7.389), thus z_AIC < z_BIC
    # In other words, AIC-correction is more prone to `model0` (complicated model); BIC-correction to `model1` (simpler model)
    
    Vuong_decision <- function(z, conf.level = .95, ddf = d_df, nm1 = nm[id1], nm0 = nm[id0]) {
      n <- length(z)
      if ((n > 1L) && (length(nm0) == 1L)) nm0 <- rep(nm0, times = n)
      if (length(nm1) != n || length(nm0) != n) stop('should not happen')
      out <- character(length = n)
      # default value: choose the simpler model
      out[idx1] <- nm1[idx1 <- (ddf < 0)]
      out[idx0] <- nm0[idx0 <- (ddf > 0)]
      out[ddf == 0] <- '[tie]'
      # Vuong test
      out[idx1] <- nm1[idx1 <- (z > qnorm(conf.level, lower.tail = TRUE))] 
      # exceeds the positive (1 − alpha)-quantile, choose model1
      out[idx0] <- nm0[idx0 <- (z < qnorm(conf.level, lower.tail = FALSE))] 
      # falls below the negative (1 − alpha)-quantile, choose model0
      return(out)
    }
    out <- data.frame(
      out0,
      z_AIC = c(NA_real_, z_AIC), z_BIC = c(NA_real_, z_BIC),
      Decision_AIC = c(NA_character_, Vuong_decision(z_AIC)),
      Decision_BIC = c(NA_character_, Vuong_decision(z_BIC)),
      check.names = FALSE
    )
    attr(out, which = 'heading') <- 'Vuong\'s closeness test'
    class(out) <- c('vuong', 'data.frame') # write ?print.vuong later
  })
  
  attr(out, which = 'logLik') <- .logL
  #attr(out, which = 'AIC') <- .AIC
  #attr(out, which = 'BIC') <- .BIC
  return(out) 
  
}

