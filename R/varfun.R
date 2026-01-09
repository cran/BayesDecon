#' Error Variance as a function of X
#'
#' @param obj A `bdeconv_result` object.
#' @param vals A numeric vector (or matrix) on which to evaluate the error
#'  variance.
#'
#' @return A numeric matrix or array of variance values.
#' @keywords internal
varfun <- function(obj, vals) {
  UseMethod('varfun')
}

#' @keywords internal
varfun.uv_res <- function(obj, vals) {
  basis <- bspline_basis(vals, obj$lwr, obj$upr, length(obj$knots))
  theta_v <- matrix(0, nrow = length(obj), ncol = length(vals))
  for (i in 1:length(obj)) {
    theta_v[i, ] <- basis %*% t(exp(obj[[i]]$theta_v)) * obj[[i]]$var_e
  }
  return(theta_v)
}

#' @keywords internal
varfun.mv_res <- function(obj, vals) {
  return(sapply(1:ncol(obj), function(dd) varfun(obj[, dd], vals[, dd])))
}
