




#' @title Sort Objects from \CRANpkg{mixsmsn} by Location Parameters
#' 
#' @description 
#' To sort an object returned from package \CRANpkg{mixsmsn} by its location parameters
#' 
#' @param x `'Normal'`, `'Skew.normal'`, `'Skew.t'` object
#' 
#' @param decreasing \link[base]{logical} scalar. Should the sort the location parameter
#' be increasing (`FALSE`, default) or decreasing (`TRUE`)?
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' \link[mixsmsn]{smsn.mix} does *not* order the location parameter
#' 
#' @returns 
#' 
#' Function [sort.Normal()] returns a `'Normal'` object.
#' 
#' Function [sort.Skew.normal()] returns a `'Skew.normal'` object.
#' 
#' Function [sort.Skew.t()] returns a `'Skew.t'` object.
#' 
#' 
#' @seealso 
#' \link[base]{sort}
#' 
#' @name sort_mixsmsn
#' @method sort Skew.normal
#' @export sort.Skew.normal
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
#' @export sort.Normal
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
#' @export sort.Skew.t
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
#' @export sort.t
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

