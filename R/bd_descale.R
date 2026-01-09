#' @keywords internal
bd_descale <- function(obj, scl, loc) {
  UseMethod('bd_descale')
}

#' @keywords internal
bd_descale.zi_res <- function(obj, scl, loc) {
  obj$lwr   <- scl * obj$lwr + loc
  obj$upr   <- scl * obj$upr + loc
  obj$knots <- scl * obj$knots + loc
  obj$w     <- scl * obj$w + loc
  obj$x     <- scl * obj$x + loc
  obj$u     <- scl * obj$u
  obj$x_nz  <- scl * obj$x_nz + loc
  obj$var_e <- obj$var_e
  obj$w_var <- scl^2 * obj$w_var
  obj$x_e <- scl * obj$x_e + loc
  return(obj)
}

#' @keywords internal
bd_descale.nc_res <- function(obj, scl, loc) {
  obj$lwr   <- scl * obj$lwr + loc
  obj$upr   <- scl * obj$upr + loc
  obj$knots <- scl * obj$knots + loc
  obj$w     <- scl * obj$w + loc
  obj$x     <- scl * obj$x + loc
  obj$u     <- scl * obj$u
  obj$x_nz  <- scl * obj$x_nz + loc
  obj$var_e <- obj$var_e
  obj$w_var <- scl^2 * obj$w_var
  obj$mu_x  <- scl * obj$mu_x + loc
  obj$sig_x <- scl * obj$sig_x + loc
  return(obj)
}

#' @keywords internal
bd_descale.mv_res <- function(obj, scl, loc) {
  obj$lwr    <- scl * obj$lwr + loc
  obj$upr    <- scl * obj$upr + loc
  obj$knots  <- scl * obj$knots + loc
  q <- obj$q
  p <- obj$p
  d <- q + p
  for (dd in 1:d) {
    obj$w[, dd, ]    <- scl[dd] * obj$w[, dd, ] + loc[dd]
    obj$x[, dd, ]    <- scl[dd] * obj$x[, dd, ] + loc[dd]
    obj$u[, dd, ]    <- scl[dd] * obj$u[, dd, ]
    obj$x_nz[, dd, ] <- scl[dd] * obj$x_nz[, dd, ] + loc[dd]
    obj$var_e[, dd]  <- obj$var_e[, dd]
    obj$w_var[dd, ]  <- scl[dd]^2 * obj$w_var[dd, ]
  }
  if (q>0){
    for (qq in 1:q) {
      obj$x_e[, qq, ] <- scl[qq] * obj$x_e[, qq, ] + loc[qq]
    }
  }
  if (p>0){
    for (pp in 1:p) {
      obj$mu_x[, pp, ]  <- scl[q + pp] * obj$mu_x[, pp, ] + loc[q + pp]
      obj$sig_x[, pp, ] <- scl[q + pp] * obj$sig_x[, pp, ] + loc[q + pp]
    }
  }
  return(obj)
}
