

#' @title Letter-Value Estimation of Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description
#' 
#' Letter-value based estimation (Hoaglin, 2006) of 
#' Tukey \eqn{g}-&-\eqn{h} distribution and its constrained versions (\eqn{g}-distribution, \eqn{h}-distribution).
#' 
#' All equation numbers mentioned below refer to Hoaglin (2006).
#' 
#' @param x \link[base]{double} vector, one-dimensional observations
#' 
#' @param p_g \link[base]{double} vector, the probabilities used for estimating parameter \eqn{g}.
#' Or, use \code{p_g = FALSE} to implement the constraint \eqn{g=0}. 
#' 
#' @param p_h \link[base]{double} vector, the probabilities used for estimating parameter \eqn{h}.
#' Or, use \code{p_h = FALSE} to implement the constraint \eqn{h=0}.
#' 
#' @param halfSpread \link[base]{character} scalar, 
#' either to use \code{'both'} half-spreads (default),
#' \code{'lower'} half-spread, or \code{'upper'} half-spread.
#' 
#' @param A,g estimated mean \eqn{\hat{A}} and skewness \eqn{\hat{g}} of Tukey's \eqn{g}-&-\eqn{h} distribution
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' \link{letterV_g} estimates parameter \eqn{g} using equation (10).
#' 
#' \link{letterV_B} estimates parameter \eqn{B} for Tukey's \eqn{g}-distribution
#' i.e., when \eqn{h=0} and \eqn{g\neq 0}, using equation (8a) and (8b).
#' 
#' \link{letterV_B_g_h} estimates parameters \eqn{B} and \eqn{h} when \eqn{g\neq 0}, using equation (33).
#' 
#' \link{letterV_B_h} estimates parameters \eqn{B} and \eqn{h} for Tukey's \eqn{h}-distribution,
#' i.e., when \eqn{g=0} and \eqn{h\neq 0}, using equation (26a), (26b) and (27).
#' 
#' \code{letterValue} plays a similar role as \code{fitdistrplus:::start.arg.default},
#' thus extends \link[fitdistrplus]{fitdist} for estimating Tukey's \eqn{g}-&-\eqn{h} distributions.
#' 
#' 
#' @return 
#' 
#' \link{letterValue} returns a \link[base]{double} vector of estimates \eqn{(\hat{A}, \hat{B}, \hat{g}, \hat{h})}
#' for a Tukey's \eqn{g}-&-\eqn{h} distribution.
#' 
#' 
#' @references 
#' Hoaglin, D.C. (2006). Summarizing Shape Numerically: The \eqn{g}-and-\eqn{h} Distributions. 
#' In \emph{Exploring Data Tables, Trends, and Shapes} (eds D.C. Hoaglin, F. Mosteller and J.W. Tukey),
#' Wiley Series in Probability and Statistics.
#' \doi{10.1002/9781118150702.ch11}
#' 
#' @seealso \link[fitdistrplus]{fitdist}
#' 
#' @examples 
#' set.seed(77652); x = rGH(n = 1e3L, g = -.3, h = .1)
#' letterValue(x, p_g = FALSE, p_h = FALSE)
#' letterValue(x, p_g = FALSE)
#' letterValue(x, p_h = FALSE)
#' 
#' (y0 = letterValue(x))
#' library(fitdistrplus)
#' fit <- fitdist(x, distr = 'GH', start = as.list.default(y0))
#' autoplot(fit)
#' 
#' 
#' @name letterValue
#' @export
letterValue <- function(
  x,
  p_g = seq.int(from = .15, to = .25, by = .005),
  p_h = seq.int(from = .15, to = .35, by = .005),
  halfSpread = c('both', 'lower', 'upper'),
  ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  raw <- list() 
  
  #### prepare parameters
  
  if (!isFALSE(p_g)) {
    if (!is.double(p_g) || anyNA(p_g)) stop('p_g (for estimating g) must be double without missing')
    p_g <- sort.int(unique_allequal(p_g))
    if (!length(p_g)) stop('`p_g` cannot be len-0')
    if (any(p_g <= 0, p_g >= .5)) stop('p_g (for estimating g) must be between 0 and .5 (not including)')
    L <- A - quantile(x, probs = p_g) # lower-half-spread (LHS), the paragraph under eq.10 on p469
    U <- quantile(x, probs = 1 - p_g) - A # upper-half-spread (UHS)
    ok <- (L != 0) & (U != 0) # may ==0 due to small sample size
    if (!all(ok)) p_g <- p_g[ok] # else do nothing
    if (!length(p_g)) stop('`p_g` cannot be len-0')
  }

  if (!isFALSE(p_h)) {
    if (!is.double(p_h) || anyNA(p_h)) stop('p_h (for estimating h) must be double without missing')
    p_h <- sort.int(unique_allequal(p_h))
    if (!length(p_h)) stop('`p_h` cannot be len-0')
    if (any(p_h <= 0, p_h >= .5)) stop('p_h (for estimating h) must be between 0 and .5 (not including)')
    L <- A - quantile(x, probs = p_h) # lower-half-spread (LHS), the paragraph under eq.10 on p469
    U <- quantile(x, probs = 1 - p_h) - A # upper-half-spread (UHS)
    ok <- (L != 0) & (U != 0) # may ==0 due to small sample size
    if (!all(ok)) p_h <- p_h[ok] # else do nothing
    if (!length(p_h)) stop('`p_h` cannot be len-0')
  }
  
  #### end of prepare parameters
  
  halfSpread <- match.arg(halfSpread)
  
  if (isFALSE(p_g)) g <- 0 else g <- raw[['g']] <- letterV_g(A = A, p_g = p_g, x = x) # must use both-spread
  
  if (isFALSE(p_h)) {
    h <- 0
    if (g == 0) return(c(A = A, B = mad(x, center = A), g = 0, h = 0)) # see ?stats::mad
    B <- raw[['B']] <- letterV_B(x, A = A, g = g, p_g = p_g, halfSpread = halfSpread)
  } else {
    Bh <- raw[['Bh']] <- if (g != 0) {
      letterV_B_g_h(A = A, g = g, p_h = p_h, x = x, halfSpread = halfSpread)
    } else letterV_B_h(A = A, p_h = p_h, x = x, halfSpread = halfSpread)
    B <- Bh$B
    h <- Bh$h
  }
  
  ret <- c(A = A, B = B, g = g, h = h)
  # attr(ret, which = 'raw') <- raw # only used by developer
  return(ret)
  
}

#' @rdname letterValue
#' @export
letterV_B_g_h <- function(A, g, p_h, x, halfSpread, ...) {
  L <- A - quantile(x, probs = p_h) # lower-half-spread (LHS), p469
  U <- quantile(x, probs = 1 - p_h) - A # upper-half-spread (UHS)
  z <- qnorm(p_h)
  regx <- z^2 / 2
  
  regyU <- log(g * U / expm1(-g*z)) # p487-eq(33)
  regyL <- log(g * L / (-expm1(g*z))) # see the equation on bottom of p486
  
  # p487, paragraph under eq(33)
  regy1 <- (regyU + regyL) / 2 # 'averaging the logrithmic results'
  regy2 <- log(g * (U+L) / (exp(-g*z) - exp(g*z))) # 'dividing the full spread with appropriate denominator', 
  # `regy2` is obtained by summing up of the numeritor and denominator of the equation on bottom of p486
  # note that (expm1(x) - expm1(y)) == (exp(x) - exp(y))
  # I prefer `regy2`
  regy <- regy2
  
  #cf1 <- lm.fit(x = cbind(1, regx), y = cbind(regy1))$coefficients
  #cf2 <- lm.fit(x = cbind(1, regx), y = cbind(regy2))$coefficients
  cf <- lm.fit(x = cbind(1, regx), y = cbind(regy))$coefficients
  cfL <- lm.fit(x = cbind(1, regx), y = cbind(regyL))$coefficients
  cfU <- lm.fit(x = cbind(1, regx), y = cbind(regyU))$coefficients
  #reg <- data.frame(x = rep(regx, times = 4L), y = c(regy1, regy2, regyL, regyU))
  reg <- data.frame(x = rep(regx, times = 3L), y = c(regy, regyL, regyU))
  
  ret <- list(
    B = exp(switch(halfSpread, both = cf[1L], upper = cfU[1L], lower = cfL[1L])),
    h = max(0, switch(halfSpread, both = cf[2L], upper = cfU[2L], lower = cfL[2L]))
    # slope `h` should be positive. When a negative slope is fitted, use h = 0
  )
  attr(ret, which = 'reg') <- reg
  attr(ret, which = 'np') <- length(p_h)
  attr(ret, which = 'A') <- A
  attr(ret, which = 'g') <- g
  #attr(ret, which = 'cf') <- cbind(both1 = cf1, both2 = cf2, lower = cfL, upper = cfU)
  attr(ret, which = 'cf') <- cbind(both = cf, lower = cfL, upper = cfU)
  class(ret) <- 'letterV_B_g_h'
  return(ret)
}


#' @rdname letterValue
#' @export
letterV_B_h <- function(A, p_h, x, halfSpread) {
  L <- A - quantile(x, probs = p_h) # lower-half-spread (LHS), p469
  U <- quantile(x, probs = 1 - p_h) - A # upper-half-spread (UHS)
  z <- qnorm(p_h)
  regx <- z^2 / 2
  
  regy <- log((U+L) / (-2*z)) # p483-eq(28); all.equal(U+L, q[-id] - q[id])
  regyL <- log(-L/z) # easy to derive from p483-eq(26a)
  regyU <- log(-U/z) # easy to derive from p483-eq(26b)
  cf <- lm.fit(x = cbind(1, regx), y = cbind(regy))$coefficients
  cfL <- lm.fit(x = cbind(1, regx), y = cbind(regyL))$coefficients
  cfU <- lm.fit(x = cbind(1, regx), y = cbind(regyU))$coefficients
  reg <- data.frame(x = rep(regx, times = 3L), y = c(regy, regyL, regyU))

  ret <- list(
    B = exp(switch(halfSpread, both = cf[1L], upper = cfU[1L], lower = cfL[1L])),
    h = max(0, switch(halfSpread, both = cf[2L], upper = cfU[2L], lower = cfL[2L]))
    # slope `h` should be positive. When a negative slope is fitted, use h = 0
  )
  attr(ret, which = 'reg') <- reg
  attr(ret, which = 'np') <- length(p_h)
  attr(ret, which = 'cf') <- cbind(both = cf, lower = cfL, upper = cfU)
  class(ret) <- 'letterV_B_h'
  return(ret)
}


#' @rdname letterValue
#' @export
letterV_B <- function(A, g, p_g, x, halfSpread) {
  # when (h == 0) && (g != 0), calculate B (p469, Example, p471-eq(11a))
  q <- quantile(x, probs = c(p_g, 1-p_g))
  z <- qnorm(p_g) # length n
  id <- seq_along(p_g)
  regx <- expm1(c(g*z, -g*z))/g # p469, eq(8a-8b)
  B <- lm.fit(x = cbind(regx), y = cbind(q - A))$coefficients
  if (B < 0) stop('estimated B < 0 ???')
  BL <- lm.fit(x = cbind(regx[id]), y = cbind(q[id] - A))$coefficients
  if (BL < 0) stop('estimated B (lower spread) < 0 ???')
  BU <- lm.fit(x = cbind(regx[-id]), y = cbind(q[-id] - A))$coefficients
  if (BU < 0) stop('estimated B (upper spread) < 0 ???')
  #labs <- c('Lower', 'Upper') # 'L' is ahead of 'U' in alphabets 
  #reg <- data.frame(x = regx, y = q - A, id = rep(labs, each = length(p_g)))
  reg <- data.frame(x = regx, y = q - A)
  ret <- unname(switch(halfSpread, both = B, upper = BU, lower = BL))
  attr(ret, which = 'B') <- c(both = B, upper = BU, lower = BL)
  attr(ret, which = 'np') <- length(p_g)
  attr(ret, which = 'reg') <- reg
  class(ret) <- 'letterV_B'
  return(ret)
}



#' @rdname letterValue
#' @export
letterV_g <- function(A, p_g, x) {
  # p469, equation (10) Estimating g; must use both half-spread
  L <- A - quantile(x, probs = p_g) # lower-half-spread (LHS), the paragraph under eq.10 on p469
  U <- quantile(x, probs = 1 - p_g) - A # upper-half-spread (UHS)
  gs <- (log(L / U) / qnorm(p_g)) # p469, equation (10)
  g <- median.default(gs) 
  # p474, 2nd paragraph, the values of $g_p$ do not show regular dependence on $p_g$,
  # thus we take their median as a constant-$g$ description of the skewness.
  attr(g, which = 'A') <- A # only for plot
  attr(g, which = 'reg') <- data.frame(gs = gs, p = p_g) # only for plot
  class(g) <- 'letterV_g'
  return(g)
}

