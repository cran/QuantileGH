


#' @title Hoaglin's Letter-Value-Based method For Estimation of 
#' Tukey \eqn{g}-&-\eqn{h} distribution and its constrained versions (\eqn{g}-distribution, \eqn{h}-distribution)
#' 
#' @description
#' 
#' \code{\link{Hoaglin_GH}} implements the letter-value-based method of estimating \code{g} and \code{h},
#' described in Hoaglin (2006), page 487, equation (33) (and the few lines beneath equation (33)).
#' 
#' @param x observations
#' 
#' @param p_g \code{FALSE} if the \eqn{g}-parameter is constrained at 0, otherwise probabilities used for estimating parameter `g`.
#' 
#' @param p_h \code{FALSE} if the \eqn{h}-parameter is constrained at 0, otherwise probabilities used for estimating parameter `h`
#' 
#' @param halfSpread \code{\link[base]{character}} scalar, which half-spread should be used
#' 
#' @param ... additional parameters, currently not in use
#' 
#' @details 
#' 
#' \code{Hoaglin_GH} can replace the functionality of \code{fitdistrplus:::start.arg.default},
#' thus extend \code{\link[fitdistrplus]{fitdist}} for estimating Tukey's \eqn{g}-&-\eqn{h} distributions.
#' 
#' All page/formula numbers refer to Hoaglin (2006).
#' 
#' @return 
#' 
#' \code{\link{Hoaglin_GH}} returns a \code{\link[base]{double}} vector with element names \code{A}, \code{B}, \code{g} and \code{h}, 
#' indicating the parameters of a Tukey's \eqn{g}-&-\eqn{h} distribution.
#' 
#' 
#' @references 
#' Hoaglin, D.C. (1006). Summarizing Shape Numerically: The \eqn{g}-and-\eqn{h} Distributions. 
#' In \emph{Exploring Data Tables, Trends, and Shapes} (eds D.C. Hoaglin, F. Mosteller and J.W. Tukey).
#' \doi{10.1002/9781118150702.ch11}
#' 
#' @examples 
#' set.seed(451); x = rGH(n = 1e2L, g = -.3, h = .1)
#' (y0 = Hoaglin_GH(x))
#' 
#' library(fitdistrplus)
#' fitdist(x, distr = 'GH', start = as.list.default(y0))
#' 
#' @export
Hoaglin_GH <- function(
  x,
  p_g = seq.int(from = .15, to = .25, by = .005),
  p_h = seq.int(from = .15, to = .35, by = .005),
  halfSpread = c('full', 'upper', 'lower', 'longer', 'shorter'),
  ...
) {
  
  if (anyNA(x)) stop('do not allow NA in observations')
  A <- median.default(x)
  halfSpread <- match.arg(halfSpread)
  # figs <- list() # no plot for end user of \pkg{QuantileGH} package
  
  if (isFALSE(p_g)) {
    Hg <- NULL
    g <- 0
  } else {
    Hg <- Hoaglin_g(x, A = A, p = p_g) # must use full-spread
    # return(Hg) # in debug
    # figs$g <- plot.Hoaglin_g(Hg)
    g <- unclass(Hg)
  }
  
  if (isFALSE(p_h)) {
    h <- 0
    if (g == 0) return(c(A = A, B = mad(x, center = A), g = 0, h = 0)) # see ?stats::mad
    B <- Hoaglin_B(x, Hg = Hg, halfSpread = halfSpread)
    #return(B) # in debug
    # figs$B <- plot.Hoaglin_B(B)
  } else {
    Bh <- Hoaglin_B_h(x, p = p_h, Hg = Hg, A = A, halfSpread = halfSpread)
    # return(Bh) # in debug
    # figs$Bh <- plot.Hoaglin_B_h(Bh)
    B <- Bh$B
    h <- Bh$h
  }
  
  # autoplot(ggstack(figs, share_legend = FALSE, nrow = 1L))
  return(c(A = A, B = unclass(B), g = g, h = h))
  
}




Hoaglin_B_h <- function(x, p, Hg, A, halfSpread) {
  if (!is.double(p) || anyNA(p) || any(p <= 0, p >= .5)) stop('p (for estimating h) must be between 0 and .5')
  if (anyDuplicated.default(p) || is.unsorted(p, strictly = TRUE)) stop('p (for estimating h) must be strictly increasing')
  has_g <- !missing(Hg) && length(Hg)
  if (has_g) A <- attr(Hg, which = 'A', exact = TRUE)
  repeat { # my addition
    id <- seq_along(p)
    q <- quantile(x, probs = c(p, 1 - p))
    L <- A - q[id] # lower-half-spread (LHS), p469 (name clash ?base::gl)
    U <- q[-id] - A # upper-half-spread (UHS)
    ok <- (L != 0) & (U != 0) # may ==0 due to small sample size
    if (all(ok)) break
    p <- p[ok]
  }
  z <- qnorm(p)
  regx <- z^2 / 2
  
  if (has_g) {
    
    if (!inherits(Hg, 'Hoaglin_g')) stop('Hg must be \'Hoaglin_g\' object')
    g <- unclass(Hg)
    if (halfSpread %in% c('longer', 'shorter')) {
      if (g < 0) {
        if (halfSpread == 'longer') halfSpread <- 'lower'
        if (halfSpread == 'shorter') halfSpread <- 'upper'
      }
      if (g > 0) {
        if (halfSpread == 'longer') halfSpread <- 'upper'
        if (halfSpread == 'shorter') halfSpread <- 'lower'
      }
    }
    regyU <- log(g * U / expm1(-g*z)) # p487-eq(33)
    regyL <- log(g * L / (-expm1(g*z))) # see the equation on bottom of p486
    # p487, paragraph under eq(33)
    regy_v1 <- (regyU + regyL) / 2 # 'averaging the logrithmic results'
    regy_v2 <- log(g * (U+L) / (exp(-g*z) - exp(g*z))) # 'dividing the full spread with appropriate denominator', 
    # `regy_v2` is obtained by summing up of the numeritor and denominator of the equation on bottom of p486
    # note that (expm1(x) - expm1(y)) == (exp(x) - exp(y))
    cf_v1 <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regy_v1))$coefficients)
    cf <- cf_v2 <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regy_v2))$coefficients)
    cfL <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regyL))$coefficients) 
    cfU <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regyU))$coefficients)
    reg <- data.frame(x = rep(regx, times = 4L), y = c(regy_v1, regy_v2, regyL, regyU), id = rep(LETTERS[1:4], each = length(p)))
    
  } else { # g == 0
    
    if (halfSpread %in% c('longer', 'shorter')) halfSpread <- 'full'
    regy <- log((U+L) / (-2*z)) # p483-eq(28); all.equal(U+L, q[-id] - q[id])
    regyL <- log(-L/z) # easy to derive from p483-eq(26a)
    regyU <- log(-U/z) # easy to derive from p483-eq(26b)
    cf <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regy))$coefficients)
    cf_v1 <- cf_v2 <- NULL
    cfL <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regyL))$coefficients)
    cfU <- as.vector(lm.fit(x = cbind(1, regx), y = cbind(regyU))$coefficients)
    reg <- data.frame(x = rep(regx, times = 3L), y = c(regy, regyL, regyU), id = rep(LETTERS[1:3], each = length(p)))
  }
  
  out <- list(
    B = exp(switch(halfSpread, full = cf[1L], upper = cfU[1L], lower = cfL[1L])),
    h = max(0, switch(halfSpread, full = cf[2L], upper = cfU[2L], lower = cfL[2L]))
    # slope `h` should be positive. When a negative slope is fitted, use h = 0
  )
  attr(out, 'p') <- p
  attr(out, 'reg') <- reg
  attr(out, 'has_g') <- has_g
  attr(out, 'cf_full') <- cf
  attr(out, 'cf_full1') <- cf_v1 
  attr(out, 'cf_full2') <- cf_v2
  attr(out, 'cf_lower') <- cfL
  attr(out, 'cf_upper') <- cfU
  class(out) <- 'Hoaglin_B_h'
  return(out)
}










Hoaglin_B <- function(x, Hg, halfSpread) {
  # when (h == 0) && (g != 0), calculate B (p469, Example, p471-eq(11a))
  if (missing(Hg) || !inherits(Hg, 'Hoaglin_g')) stop('Hg must be \'Hoaglin_g\' object')
  Hg_attr <- attributes(Hg)
  g <- unclass(Hg)
  A <- Hg_attr[['A']]
  q <- Hg_attr[['q']] # length 2n
  p <- Hg_attr[['p']] # length n
  z <- qnorm(p) # length n
  id <- seq_along(p)
  regx <- expm1(c(g*z, -g*z))/g # p469, eq(8a-8b)
  B <- as.vector(lm.fit(x = cbind(regx), y = cbind(q - A))$coefficients)
  if (B < 0) stop('estimated B < 0 ???')
  BL <- as.vector(lm.fit(x = cbind(regx[id]), y = cbind(q[id] - A))$coefficients)
  if (BL < 0) stop('estimated B (lower spread) < 0 ???')
  BU <- as.vector(lm.fit(x = cbind(regx[-id]), y = cbind(q[-id] - A))$coefficients)
  if (BU < 0) stop('estimated B (upper spread) < 0 ???')
  labs <- c('Lower', 'Upper') # 'L' is ahead of 'U' in alphabets 
  reg <- data.frame(x = regx, y = q - A, id = rep(labs, each = length(p)))
  if (halfSpread %in% c('longer', 'shorter')) {
    if (g < 0) {
      if (halfSpread == 'longer') halfSpread <- 'lower'
      if (halfSpread == 'shorter') halfSpread <- 'upper'
    }
    if (g > 0) {
      if (halfSpread == 'longer') halfSpread <- 'upper'
      if (halfSpread == 'shorter') halfSpread <- 'lower'
    }
  }
  out <- switch(halfSpread, full = B, upper = BU, lower = BL)
  attr(out, 'B_full') <- B 
  attr(out, 'B_upper') <- BU
  attr(out, 'B_lower') <- BL
  attr(out, 'reg') <- reg
  class(out) <- 'Hoaglin_B'
  return(out)
}








Hoaglin_g <- function(x, A = median.default(x), p = seq.int(from = 0, to = .45, by = .005)) {
  # p468, Estimating g; must use both half-spread
  if (!is.double(p) || anyNA(p)) stop('p (for estimating g) must be double without missing')
  if (any(p <= 0, p >= .5)) stop('p (for estimating g) must be between 0 and .5 (not including)')
  if (anyDuplicated.default(p) || is.unsorted(p, strictly = TRUE)) stop('p (for estimating g) must be strictly increasing')
  repeat { # mine addition
    id <- seq_along(p)
    q <- quantile(x, probs = c(p, 1 - p))
    L <- A - q[id] # lower-half-spread (LHS), p469
    U <- q[-id] - A # upper-half-spread (UHS)
    ok <- (L != 0) & (U != 0) # may ==0 due to small sample size
    if (all(ok)) break
    p <- p[ok]
  }
  gs <- (log(L / U) / qnorm(p)) # $g_p$ of p469-eq(10) (also p487-eq(31))
  ghat <- median.default(gs) # is this the best choice?
  attr(ghat, 'A') <- A
  attr(ghat, 'gs') <- gs
  attr(ghat, 'p') <- p
  attr(ghat, 'q') <- q
  class(ghat) <- 'Hoaglin_g'
  return(ghat)
}

