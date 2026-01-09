#' Density deconvolution using bayesian semiparametric methods
#' @param x A `bdeconv_result` object.
#' @param ... Unused.
#' @keywords internal
as.list.bdeconv_result <- function(x, ...) {
  return(as.list(unclass(x)))
}
