

#' @title Quantile Least Mahalanobis Distance Estimator for Tukey \eqn{g}-&-\eqn{h} Mixture
#'
#' @description
#'
#' \tabular{ll}{
#' Package: \tab QuantileGH\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1.0\cr
#' Date: \tab 2021-12-15\cr
#' License: \tab GPL (>= 2)\cr
#' }
#'
#' @details 
#' 
#' This package provides tools for simulating and fitting finite mixtures of the 4-parameter Tukey \eqn{g}-&-\eqn{h} distributions. 
#' Tukey \eqn{g}-&-\eqn{h} mixture is highly flexible to model multimodal distributions with variable degree of skewness and kurtosis in the components. 
#' The Quantile Least Mahalanobis Distance estimator (\code{\link{QLMDe}}) is used for estimating parameters of the finite Tukey \eqn{g}-&-\eqn{h} mixtures.
#' \code{\link{QLMDe}} is an indirect estimator that minimizes the Mahalanobis distance between the sample and model-based quantiles.
#' A backward-forward stepwise model selection algorithm is provided to find
#' \itemize{
#' \item {a parsimonious Tukey \eqn{g}-&-\eqn{h} mixture model, conditional on a given number-of-components; and}
#' \item {the optimal number of components within the user-specified range.}
#' }
#'
# @references none yet
#'
#'
#' @examples
#' # see ?QLMDe
#'
#' @import graphics methods stats utils
#' @import ggplot2
#' 
#' @importFrom goftest cvm.test
#' @importFrom LaplacesDemon KLD
#' @importFrom rstpm2 vuniroot
#' @importFrom scales percent
#' @importFrom tclust tkmeans
#'
#' @docType package
#' @keywords package
#' @name QuantileGH-package
#' @aliases QuantileGH
NULL
