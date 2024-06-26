% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step_fmx.R
\name{step_fmx}
\alias{step_fmx}
\title{Forward Selection of \eqn{gh}-parsimonious Model with Fixed Number of Components \eqn{K}}
\usage{
step_fmx(
  object,
  test = c("BIC", "AIC"),
  direction = c("forward", "backward"),
  ...
)
}
\arguments{
\item{object}{\linkS4class{fmx} object}

\item{test}{\link[base]{character} scalar, criterion to be used, either
Akaike's information criterion \link[stats]{AIC}-like, or
Bayesian information criterion \link[stats]{BIC}-like (default).}

\item{direction}{\link[base]{character} scalar, \code{'forward'} (default) or
\code{'backward'}}

\item{...}{additional parameters, currently not in use}
}
\value{
Function \link{step_fmx} returns an object of S3 class \code{'step_fmx'},
which is a \link[base]{list} of selected models (in reversed order) with attribute(s)
\code{'direction'} and
\code{'test'}.
}
\description{
To select the \eqn{gh}-parsimonious mixture model,
i.e., with some \eqn{g} and/or \eqn{h} parameters equal to zero,
conditionally on a fixed number of components \eqn{K}.
}
\details{
The algorithm starts with quantile least Mahalanobis distance estimates
of either the full mixture of Tukey \eqn{g}-&-\eqn{h} distributions model, or
a constrained model (i.e., some \eqn{g} and/or \eqn{h} parameters equal to zero according to the user input).
Next, each of the non-zero \eqn{g} and/or \eqn{h} parameters is tested using the likelihood ratio test.
If all tested \eqn{g} and/or \eqn{h} parameters are significantly different from zero at the level 0.05
the algorithm is stopped and the initial model is considered \eqn{gh}-parsimonious.
Otherwise, the \eqn{g} or \eqn{h} parameter with the largest p-value is constrained to zero
for the next iteration of the algorithm.

The algorithm iterates until only significantly-different-from-zero \eqn{g} and \eqn{h} parameters
are retained, which corresponds to \eqn{gh}-parsimonious Tukey \eqn{g}-&-\eqn{h} mixture model.
}
\seealso{
\link[stats]{step}
}
