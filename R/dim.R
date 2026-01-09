#' @exportS3Method
dim.bdeconv_result <- function(x) {
  obj <- x
  return(c(length(obj), obj$p + obj$q))
}
