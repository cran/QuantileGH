
# use 

#sn::dsn
#sn::psn
#sn::qsn
#sn::rsn

# and 

#sn::dst
#sn::pst
#sn::qst
#sn::rsn

# instead

#' @title Skew Normal as in \CRANpkg{mixsmsn}
#' 
#' @description ..
#' 
#' @param x,q ..
#' 
#' @param mean ..
#' 
#' @param sd ..
#' 
#' @param shape ..
#' 
#' @param lower.tail ..
#' 
#' @param log,log.p ..
#' 
#' @seealso \code{mixsmsn:::dSN} \link[sn]{dsn}
#' 
#' @return 
#' \link{dSN} returns a \link[base]{numeric} \link[base]{vector} of the same length as \code{x},
#' or a \link[base]{numeric} \link[base]{matrix} of the same dimension as \code{x}.
#' 
#' @name mixsmsn_skewNormal
#' @export
dSN <- function(x, mean = 0, sd = 1, shape = 0, log = FALSE) {
  logd <- log(2) + dnorm(x, mean = mean, sd = sd, log = TRUE) + pnorm(shape*(x-mean)/sd, log.p = TRUE)
  if (log) return(logd)
  return(exp(logd))
}

#' @rdname mixsmsn_skewNormal
#' @export
pSN <- function(q, mean = 0, sd = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  stop('not written yet')
}


#' @title Skew \eqn{t} as in \CRANpkg{mixsmsn}
#' 
#' @description ..
#' 
#' @param x,q ..
#' 
#' @param mean ..
#' 
#' @param sd ..
#' 
#' @param shape ..
#' 
#' @param nu ..
#' 
#' @param lower.tail ..
#' 
#' @param log,log.p ..
#' 
#' @seealso \code{mixsmsn:::dt.ls} \link[sn]{dst}
#' 
#' @return
#' 
#' \link{dST} returns a \link[base]{numeric} \link[base]{vector} of the same length as \code{x},
#' or a \link[base]{numeric} \link[base]{matrix} of the same dimension as \code{x}.
#' 
#' @name mixsmsn_skewT
#' @export
dST <- function(x, mean = 0, sd = 1, shape = 0, nu = 4, log = FALSE) {
  e <- (x - mean)/sd
  logd <- log(2) + dt(e, df = nu, log = TRUE) + pt(q = sqrt((1 + nu)/(e^2 + nu)) * e * shape, df = 1 + nu, log.p = TRUE) - log(sd)
  if (log) return(logd)
  return(exp(logd))
}


#' @name mixsmsn_skewT
#' @export
pST <- function(q, mean = 0, sd = 1, shape = 0, nu = 4, lower.tail = TRUE, log.p = FALSE) {
  stop('not written yet')
}