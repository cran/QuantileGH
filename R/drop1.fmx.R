


#' @title Drop or Add One Parameter from \linkS4class{fmx} Object
#' 
#' @description 
#' 
#' Fit \linkS4class{fmx} models with a single parameters being added or dropped.
#' 
#' @param object \linkS4class{fmx} object
#' 
#' @param ... additional parameters, currently not in use.
#' 
#' @details ..
#' 
#' @note 
#' Note that \link{drop1.fmx} and \link{add1.fmx} do \strong{not} return an \link[stats]{anova} table, like other
#' \code{stats:::drop.*} or \code{stats:::add1.*} functions do.
#' 
#' @return
#' 
#' \link{drop1.fmx} and \link{add1.fmx} return a \link[base]{list} of \linkS4class{fmx} objects,
#' in the reverse order of model selection.
#' 
#' @seealso \link[stats]{step} \link[stats]{drop1} \link[stats]{add1}
#' 
#' @importFrom stats drop1 add1
#' 
#' @examples 
#' 
#' \donttest{ 
#' 
#' # donttest to save time
#' 
#' (d2 = fmx('GH', A = c(1,6), B = 1.2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
#' set.seed(3123); hist(x2 <- rfmx(n = 1e3L, dist = d2))
#' system.time(m0 <- QLMDe(x2, distname = 'GH', K = 2L, constraint = c('g1', 'g2', 'h1', 'h2')))
#' system.time(m1 <- QLMDe(x2, distname = 'GH', K = 2L, constraint = c('g1', 'h2')))
#' system.time(m2 <- QLMDe(x2, distname = 'GH', K = 2L)) # ~2 secs
#' 
#' d1 = drop1(m1)
#' d1 # NULL
#' d2 = drop1(m2)
#' vapply(d2, FUN = constraint_TeX, FUN.VALUE = '')
#' 
#' a0 = add1(m0)
#' vapply(a0, FUN = constraint_TeX, FUN.VALUE = '')
#' a1 = add1(m1)
#' vapply(a1, FUN = constraint_TeX, FUN.VALUE = '')
#' 
#' }
#' 
#' @name drop1_fmx
#' @export
drop1.fmx <- function(object, ...) {
  if (!length(object@data)) return(invisible())
  K <- dim(object@pars)[1L]
  
  # existing constraint of `object`, to be respected
  constr_ <- attr(fmx_constraint(object), which = 'user', exact = TRUE)
  
  # candidate parameter(s) to be dropped
  candpar <- setdiff(x = switch(object@distname, GH = {
    c(outer(c('g', 'h'), 1:K, FUN = paste0))
  }, character()), y = constr_)
  if (!length(candpar)) return(invisible()) # exception handling
  
  mods0 <- lapply(candpar, FUN = function(i) {
    i <- c(constr_, i)
    message(paste0(object@distname, K), ' ', paste(c(i, '0'), collapse = '='), ' .. ', appendLF = FALSE)
    ret <- QLMDe(object@data, data.name = object@data.name, distname = object@distname, K = K, probs = object@probs, constraint = i)
    message('done!')
    return(ret)
  })
  return(mods0[lengths(mods0) > 0L]) # in case any drop1-model runs into error
  # len-0 list is ok
}






# need to re-use the examples, to save time!
#' @rdname drop1_fmx
#' @export
add1.fmx <- function(object, ...) {
  if (!length(object@data)) return(invisible())
  K <- dim(object@pars)[1L]
  
  # existing constraint of `object` (i.e., candidate parameter(s) to be added)
  constr_ <- attr(fmx_constraint(object), which = 'user', exact = TRUE)
  if (!length(constr_)) return(invisible()) # exception handling
  
  mods0 <- lapply(seq_along(constr_), FUN = function(i) {
    i <- constr_[-i] # remove one constraint parameter
    message(paste0(object@distname, K), ' ', if (length(i)) paste(c(i, '0'), collapse = '='), ' .. ', appendLF = FALSE)
    ret <- QLMDe(object@data, data.name = object@data.name, distname = object@distname, K = K, probs = object@probs, constraint = i)
    message('done!')
    return(ret)
  })
  return(mods0[lengths(mods0) > 0L]) # in case any drop1-model runs into error
  # len-0 list is ok
}


