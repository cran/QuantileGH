




#' @title Sort Objects from \pkg{mixsmsn} by Location Parameters
#' 
#' @description 
#' To sort an object returned from package \pkg{mixsmsn} by its location parameters
#' 
#' @param x \code{Normal}, \code{Skew.normal}, \code{Skew.t} object
#' 
#' @param decreasing \link[base]{logical} scalar. Should the sort the location parameter
#' be increasing (\code{FALSE}, default) or decreasing (\code{TRUE})?
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' \link[mixsmsn]{smsn.mix} does \strong{not} order the location parameter
#' 
#' @return 
#' 
#' \link{sort.Normal} returns a \code{Normal} object.
#' 
#' \link{sort.Skew.normal} returns a \code{Skew.normal} object.
#' 
#' \link{sort.Skew.t} returns a \code{Skew.t} object.
#' 
#' 
#' @seealso \link[base]{sort}
#' 
#' @name sort_mixsmsn
#' @method sort Skew.normal
#' @export
sort.Skew.normal <- function(x, decreasing = FALSE, ...) {
  # stop on multivariable object ..
  ret <- x
  loc <- x[['mu']]
  
  o <- order(loc, decreasing = decreasing)
  ret[['mu']] <- x[['mu']][o]
  ret[['sigma2']] <- x[['sigma2']][o]
  ret[['shape']] <- x[['shape']][o]
  ret[['pii']] <- x[['pii']][o]
  if (length(x[['group']])) ret[['group']] <- match(x[['group']], table = o) # wow!
  return(ret)
}



#' @rdname sort_mixsmsn
#' @export
sort.Normal <- function(x, decreasing = FALSE, ...) {
  # stop on multivariable object ..
  ret <- x
  loc <- x[['mu']]
  
  o <- order(loc, decreasing = decreasing)
  ret[['mu']] <- x[['mu']][o]
  ret[['sigma2']] <- x[['sigma2']][o]
  ret[['pii']] <- x[['pii']][o]
  if (length(x[['group']])) ret[['group']] <- match(x[['group']], table = o) # wow!
  return(ret)
}



#' @rdname sort_mixsmsn
#' @method sort Skew.t
#' @export
sort.Skew.t <- function(x, decreasing = FALSE, ...) {
  # stop on multivariable object ..
  ret <- x
  loc <- x[['mu']]
  
  o <- order(loc, decreasing = decreasing)
  ret[['mu']] <- x[['mu']][o]
  ret[['sigma2']] <- x[['sigma2']][o]
  ret[['shape']] <- x[['shape']][o]
  if (length(x[['nu']]) != 1L) stop('mixsmsn package update?')
  ret[['pii']] <- x[['pii']][o]
  if (length(x[['group']])) ret[['group']] <- match(x[['group']], table = o) # wow!
  return(ret)
}



#' @rdname sort_mixsmsn
#' @export
sort.t <- function(x, decreasing = FALSE, ...) {
  # stop on multivariable object ..
  ret <- x
  loc <- x[['mu']]
  
  o <- order(loc, decreasing = decreasing)
  ret[['mu']] <- x[['mu']][o]
  ret[['sigma2']] <- x[['sigma2']][o]
  if (length(x[['nu']]) != 1L) stop('mixsmsn package update?')
  ret[['pii']] <- x[['pii']][o]
  if (length(x[['group']])) ret[['group']] <- match(x[['group']], table = o) # wow!
  return(ret)
}







#' @export
logLik.Skew.normal <- function(object, data = stop('must provide data explicitly'), ...) {
  logLik.fmx(as.fmx.Skew.normal(object, data = data, ...))
}

#' @export
logLik.Skew.t <- function(object, data = stop('must provide data explicitly'), ...) {
  logLik.fmx(as.fmx.Skew.t(object, data = data, ...))
}

