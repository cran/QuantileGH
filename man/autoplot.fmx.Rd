% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_fmx.R
\name{autoplot.fmx}
\alias{autoplot.fmx}
\title{Plot \linkS4class{fmx} Objects using \CRANpkg{ggplot2}}
\usage{
\method{autoplot}{fmx}(
  object,
  xlab = attr(object, which = "data.name", exact = TRUE),
  ylab = NULL,
  title = TeX(getTeX(object)),
  caption = NULL,
  ...
)
}
\arguments{
\item{object}{\linkS4class{fmx} object}

\item{xlab, ylab, title, caption}{\link[base]{character} scalars, the
horizontal and vertical label, title and caption}

\item{...}{potential parameters of \link{autolayer_fmx_continuous} and \link{autolayer_fmx_discrete}}
}
\value{
\link{autoplot.fmx} returns a \link[ggplot2]{ggplot} object.
}
\description{
Plot \linkS4class{fmx} objects using \CRANpkg{ggplot2}.
}
\examples{
(d2 = fmx('GH', A = c(1,6), B = 2, g = c(0,.3), h = c(.2,0), w = c(1,2)))
curve(dfmx(x, dist = d2), xlim = c(-3, 11))
curve(pfmx(x, dist = d2), xlim = c(-3, 11))
autoplot(d2)
autoplot(d2, type = 'distribution')

}
\seealso{
\link{autolayer_fmx_continuous} \link{autolayer_fmx_discrete}
\link[ggplot2]{autoplot}
}
