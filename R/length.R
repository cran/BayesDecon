#' @exportS3Method
length.bdeconv_result <- function(x) {
  if (is.null(dim(x$x))) {
    return(1)
  } else if (length(dim(x$x)) == 2) {
    if ((x$q + x$p) > 1) {
      return(1)
    }
  }
  return(dim(x$x)[[1]])
}

