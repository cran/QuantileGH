


vuniroot2 <- function(
  y, yok = !is.na(y), f, interval = stop('must provide interval'), 
  else_return = stop(err_message),
  ...
) {
  
  if (any(is.infinite(y))) stop('infinite return from function `f` wont happen') # may mask later
  out <- rep(NA_real_, times = length(y)) # y * NA_real_; # since NA * Inf -> NaN
  out[yok] <- Inf
  yok_ <- y[yok]
  
  lower <- interval[1L]
  upper <- interval[2L]
  f.interval <- f(interval) # len-2
  if (f.interval[1L] == f.interval[2L]) {
    err_message <- 'function in `vuniroot2` not monotone'
    return(else_return)
  }
  if (anyNA(f.interval)) {
    #f <<- f; interval <<- interval
    cat('Interval', paste(interval, collaspe = ', '), '\n')
    cat(f(interval), '\n')
    stop('parametrization has been optimized! Let `vuniroot2` err.')
  }
  
  #if (any(is.infinite(f.interval))) {
  #  # \code{fitdistrplus:::test1fun} requires NaN output
  #  out[yok] <- NaN 
  #  return(out)
  #}
  
  f.lower <- f.interval[1L] - yok_
  f.upper <- f.interval[2L] - yok_
  lower_neg <- (f.lower < 0)
  upper_pos <- (f.upper > 0)
  
  id_same_sign <- xor(lower_neg, upper_pos) # '==' case considered :))
  
  if (any(id_same_sign)) {
    out[yok][id_same_sign & (abs(f.lower) < abs(f.upper))] <- -Inf
    id_oppo_sign <- which(!id_same_sign) # remove .Internal for CRAN
  } else {
    id_oppo_sign <- seq_along(yok_)
  }
  
  nn <- length(id_oppo_sign)
  suppressWarnings(
    out[yok][id_oppo_sign] <- vuniroot(
      f = \(x) f(x) - yok_[id_oppo_sign],
      lower = rep(lower, times = nn), upper = rep(upper, times = nn), 
      f.lower = f.lower[id_oppo_sign], f.upper = f.upper[id_oppo_sign], ...)[[1L]]
  )
  
  return(out)
  
}

