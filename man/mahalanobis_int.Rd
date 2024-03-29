% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mahalanobis_int.R
\name{mahalanobis_int}
\alias{mahalanobis_int}
\title{A Simpler and Faster Mahalanobis Distance}
\usage{
mahalanobis_int(x, center, invcov)
}
\arguments{
\item{x}{\link[base]{numeric} \link[base]{vector}}

\item{center}{\link[base]{numeric} \link[base]{vector}, mean \eqn{\mathbf{\mu}}}

\item{invcov}{\link[base]{numeric} \link[base]{matrix}, \emph{inverted} variance-covariance \eqn{\mathbf{\Sigma}}}
}
\value{
\link{mahalanobis_int} returns a \link[base]{numeric} scalar.
}
\description{
A simpler and faster \link[stats]{mahalanobis} distance.
}
