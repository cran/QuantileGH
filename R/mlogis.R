
#' @title Transformation between Multinomial Probabilities & Logits
#' 
#' @description 
#' 
#' Performs transformation between vectors of multinomial probabilities and multinomial logits.
#' 
#' This transformation is a generalization of \code{\link[stats]{plogis}} which converts scalar logit into probability
#' and \code{\link[stats]{qlogis}} which converts probability into scalar logit.
#' 
#' @param p 'numeric' vector of multinomial probabilities, adding up to 1
#' 
#' @param q 'numeric' vector of multinomial logits
#' 
#' @details 
#' 
#' \code{\link{pmlogis_first}} and \code{\link{pmlogis_last}} take a length \eqn{k-1} 'numeric' vector of 
#' multinomial logits and convert them into length \eqn{k} multinomial probabilities, regarding the first or last category 
#' as reference, respectively.
#' 
#' \code{\link{qmlogis_first}} and \code{\link{qmlogis_last}} take a length \eqn{k} 'numeric' vector of 
#' multinomial probabilities and convert them into length \eqn{k-1} multinomial logits, regarding the first or last category 
#' as reference, respectively.
#' 
#' @return 
#' 
#' \code{\link{pmlogis_first}} and \code{\link{pmlogis_last}} return a vector of multinomial probabilities.
#' 
#' \code{\link{qmlogis_first}} and \code{\link{qmlogis_last}} returns a vector of multinomial logits, regarding the first or last category 
#' as reference, respectively.
#'   
#' @examples
#' (a = qmlogis_last(c(2,5,3)))
#' (b = qmlogis_first(c(2,5,3)))
#' pmlogis_last(a)
#' pmlogis_first(b)
#' 
#' q0 = .8300964
#' (p1 = pmlogis_last(q0))
#' (q1 = qmlogis_last(p1))
#' 
#' # various exceptions
#' pmlogis_first(qmlogis_first(c(1, 0)))
#' pmlogis_first(qmlogis_first(c(0, 1)))
#' pmlogis_first(qmlogis_first(c(0, 0, 1)))
#' pmlogis_first(qmlogis_first(c(0, 1, 0, 0)))
#' pmlogis_first(qmlogis_first(c(1, 0, 0, 0)))
#' pmlogis_last(qmlogis_last(c(1, 0)))
#' pmlogis_last(qmlogis_last(c(0, 1)))
#' pmlogis_last(qmlogis_last(c(0, 0, 1)))
#' pmlogis_last(qmlogis_last(c(0, 1, 0, 0)))
#' pmlogis_last(qmlogis_last(c(1, 0, 0, 0)))
#' 
#' @name mlogis
#' @export
qmlogis_first <- function(p) {
  lp <- log(p)
  lp[-1L] - lp[1L]
}

#' @rdname mlogis
#' @export
qmlogis_last <- function(p) {
  lp <- log(p)
  n <- length(p)
  lp[-n] - lp[n]
}



#' @rdname mlogis
#' @export
pmlogis_first <- function(q) {
  eq <- exp(q)
  id <- is.infinite(eq)
  if (any(id)) {
    if (all(id) && all(q < 0)) return(c(1, q*0))
    if (sum(id) > 1L || q[id] < 0) stop('must be single positive Inf')
    return(c(0, id))
  }
  y0 <- c(1, eq)
  return(y0/sum(y0))
}

#' @rdname mlogis
#' @export
pmlogis_last <- function(q) {
  eq <- exp(q)
  id <- is.infinite(eq)
  if (any(id)) {
    if (all(id) && all(q < 0)) return(c(q*0, 1))
    if (sum(id) > 1L || q[id] < 0) stop('must be single positive Inf')
    return(c(id, 0))
  }
  y0 <- c(eq, 1)
  return(y0/sum(y0))
}



