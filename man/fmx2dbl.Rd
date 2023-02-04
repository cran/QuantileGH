% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmx_repar.R
\name{fmx2dbl}
\alias{fmx2dbl}
\title{Reparameterization of \linkS4class{fmx} Object}
\usage{
fmx2dbl(
  x,
  distname = x@distname,
  pars = x@pars,
  K = dim(pars)[1L],
  w = x@w,
  ...
)
}
\arguments{
\item{x}{\linkS4class{fmx} object}

\item{distname}{\link[base]{character} scalar, default value from \code{x@distname}}

\item{pars}{\link[base]{numeric} \link[base]{matrix}, default value from \code{x@pars}}

\item{K}{\link[base]{integer} scalar, default value from \code{x}}

\item{w}{\link[base]{numeric} \link[base]{vector}, default value from \code{x@w}}

\item{...}{additional parameters, not currently used}
}
\value{
\link{fmx2dbl} returns a \link[base]{numeric} vector
}
\description{
To convert the parameters of \linkS4class{fmx} object into unrestricted parameters.
}
\details{
For the first parameter
\itemize{
\item {\eqn{A_1 \rightarrow A_1}}
\item {\eqn{A_2 \rightarrow A_1 + \exp(\log(d_1))}}
\item {\eqn{A_k \rightarrow A_1 + \exp(\log(d_1)) + \cdots + \exp(\log(d_{k-1}))}}
}

For mixing proportions to multinomial logits.

For \code{'norm'}: \code{sd -> log(sd)}
for \code{'GH'}: \code{B -> log(B), h -> log(h)}
}
\seealso{
\link{dbl2fmx}
}