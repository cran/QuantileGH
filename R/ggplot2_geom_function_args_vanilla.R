
#' @title \link[ggplot2]{geom_function} with Multiple Sets of Arguments
#' 
#' @description ..
#' 
#' @param args *named* \link[base]{list} of arguments
#' 
#' @param ... parameters of \link[ggplot2]{geom_function}, 
#' most importantly the parameter `fun`
#' 
#' @details 
#' Function [geom_function_args()] plots *one* function using a *list of* arguments, 
#' by calling \link[ggplot2]{geom_function} repetitively.
#' The colour labels are the names of argument list `args`.
#' 
#' @returns 
#' Function [geom_function_args()] returns a \link[base]{list} of \link[ggplot2]{ggplot} \link[ggplot2]{layer}s.
#' 
#' @note
#' Parameter `args` of \link[ggplot2]{geom_function} is *not* vectorized.
#' 
#' See \link[ggplot2]{geom_function}, for the difference from \link[ggplot2]{stat_function}.
#' 
#' @examples
#' ggplot() + 
#'  geom_function_args(
#'    args = c('$\\alpha$' = 1, '$\\beta$' = 2), 
#'    fun = function(x, a) a*x^2, 
#'    xlim = c(-3, 3)) + 
#'  labs(colour = 'Args')
#' 
#' @importFrom ggplot2 geom_function scale_colour_discrete
#' @importFrom latex2exp TeX
#' @export
geom_function_args <- function(args, ...) {
  nms <- names(args)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) stop('args must be fully named')
  if ((nag <- length(nms)) > 99L) stop('only plot up to 99 curves on one figure')
  aseq <- seq_len(nag)
  achr <- sprintf(fmt = '%02d', aseq)
  ret <- lapply(aseq, FUN = function(i) geom_function(mapping = aes(colour = achr[i]), args = args[[i]], ...))
  return(c(ret, list(
    scale_colour_discrete(breaks = achr, labels = unname(TeX(nms))) # restore `nms`
    # `unname()` is necessary!  see
    # https://stackoverflow.com/questions/44310088/how-to-add-latex-code-in-ggplot2-legend-labels
  )))
}