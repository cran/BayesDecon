#' @export
`[[.bdeconv_result` <- function(obj, i, j, ..., exact = TRUE) {
  if (!exact) {
    warning('`exact` is unimplemented for `bdeconv_result`')
  }
  if (length(i) < 1) {
    stop('attempt to select less than one element in get1index')
  } else if (length(i) > 1) {
    stop('subscript out of bounds')
  }
  if (i > length(obj)) {
    stop('subscript out of bounds')
  }
  if (length(obj) != 1) {
    if (dim(obj)[[2]] == 1) {
      if (!missing(j)) {
        stop('incorrect number of dimensions')
      }
      obj <- obj[i]
    } else {
      obj <- obj[i, j, ..., drop = FALSE]
    }
  }
  return(obj)
}


#' Subset Methods for `bdeconv_result` object
#' @param obj A `bdeconv_result` object.
#' @param i An integer row index (or indices) to select some specific MCMC iterations.
#' @param j And integer column index (or indices) to select some specific set of variables (For multivariate case only).
#' @name subset_result
#' @rdname subset_result
#' @return A `bdeconv_result` object containing the selected iterations of posterior samples for the specified variables (in the multivariate case), suitable for further analysis and visualization.
#' @export
`[.zi_res` <- function(obj, i) {
  obj$w <- obj$w[i, , drop = FALSE]
  obj$x <- obj$x[i, , drop = FALSE]
  obj$theta_x <- obj$theta_x[i, , drop = FALSE]
  obj$x_e <- obj$x_e[i, , drop = FALSE]
  obj$var_e <- obj$var_e[i]
  obj$theta_e <- obj$theta_e[i, , drop = FALSE]
  obj$x_nz <- obj$x_nz[i, , drop = FALSE]
  obj$theta_v <- obj$theta_v[i, , drop = FALSE]
  obj$theta_dens <- obj$theta_dens[i, , drop = FALSE]
  obj$u <- obj$u[i, , drop = FALSE]
  obj$k_u <- obj$k_u[i]
  obj$z_u <- obj$z_u[i, , drop = FALSE]
  obj$pi_u <- obj$pi_u[i, , drop = FALSE]
  obj$p_u <- obj$p_u[i, , drop = FALSE]
  obj$mu_u <- obj$mu_u[i, , drop = FALSE]
  obj$sig1_u <- obj$sig1_u[i, , drop = FALSE]
  obj$sig2_u <- obj$sig2_u[i, , drop = FALSE]
  return(obj)
}

#' @rdname subset_result
#' @export
`[.nc_res` <- function(obj, i) {
  obj$x <- obj$x[i, , drop = FALSE]
  obj$z_x <- obj$z_x[i, , drop = FALSE]
  obj$pi_x <- obj$pi_x[i, , drop = FALSE]
  obj$mu_x <- obj$mu_x[i, , drop = FALSE]
  obj$sig_x <- obj$sig_x[i, , drop = FALSE]
  obj$theta_v <- obj$theta_v[i, , drop = FALSE]
  obj$var_e <- obj$var_e[i]
  obj$u <- obj$u[i, , drop = FALSE]
  obj$k_u <- obj$k_u[i]
  obj$z_u <- obj$z_u[i, , drop = FALSE]
  obj$pi_u <- obj$pi_u[i, , drop = FALSE]
  obj$p_u <- obj$p_u[i, , drop = FALSE]
  obj$mu_u <- obj$mu_u[i, , drop = FALSE]
  obj$sig1_u <- obj$sig1_u[i, , drop = FALSE]
  obj$sig2_u <- obj$sig2_u[i, , drop = FALSE]
  return(obj)
}

#' @rdname subset_result
#' @export
`[.mv_res` <- function(obj, i, j) {
  # Subset i
  if (missing(i)) {
    i <- 1:length(obj)
  } else {
    obj$w <-       obj$w[i, , , drop = FALSE]
    obj$x <-       obj$x[i, , , drop = FALSE]
    if (!is.null(obj$z_x)) {
      obj$z_x <- obj$z_x[i, , , drop = FALSE]
    }
    obj$pi_x <-    obj$pi_x[i, , , drop = FALSE]
    obj$mu_x <-    obj$mu_x[i, , , drop = FALSE]
    obj$sig_x <-   obj$sig_x[i, , , drop = FALSE]
    obj$theta_x <- obj$theta_x[i, , , drop = FALSE]
    obj$theta_v <- obj$theta_v[i, , , drop = FALSE]
    obj$x_nz <-    obj$x_nz[i, , , drop = FALSE]
    obj$x_e <-     obj$x_e[i, , , drop = FALSE]
    obj$theta_e <- obj$theta_e[i, , , drop = FALSE]
    obj$var_e <-   obj$var_e[i, , drop = FALSE]
    obj$u <-       obj$u[i, , , drop = FALSE]
    obj$k_u <-     obj$k_u[i, , drop = FALSE]
    if (!is.null(obj$z_u)) {
      obj$z_u <- obj$z_u[i, , , drop = FALSE]
    }
    obj$pi_u <-    obj$pi_u[i, , , drop = FALSE]
    obj$p_u <-     obj$p_u[i, , drop = FALSE]
    obj$mu_u <-    obj$mu_u[i, , drop = FALSE]
    obj$sig1_u <-  obj$sig1_u[i, , drop = FALSE]
    obj$sig2_u <-  obj$sig2_u[i, , drop = FALSE]
    obj$corr_x <-  obj$corr_x[i, , , drop = FALSE]
    obj$corr_e <-  obj$corr_e[i, , , drop = FALSE]
  }
  # Subset j
  if (missing(j)) {
    j <- 1:(obj$q + obj$p)
    if (obj$q != 0) {
      jq <- 1:obj$q
    } else {
      jq <- integer(0)
    }
    if (obj$p != 0) {
      jp <- 1:obj$p
    } else {
      jp <- integer(0)
    }
  } else {
    jq <- j[j <= obj$q]
    jp <- j[j > obj$q] - obj$q
  }
  if (length(j) == 1) {
    if (length(jq) == 1) {
      obj$w <-       obj$w      [, jq, , drop = TRUE]
      obj$x <-       obj$x      [, jq, , drop = TRUE]
      obj$theta_x <- obj$theta_x[, jq, , drop = TRUE]
      obj$theta_v <- obj$theta_v[, jq, , drop = TRUE]
      obj$theta_e <- obj$theta_e[, jq, , drop = TRUE]
      obj$var_e <-   obj$var_e  [, jq,   drop = TRUE]
      obj$x_e <-     obj$x_e    [, jq, , drop = TRUE]
      obj$x_nz <-    obj$x_nz   [, jq, , drop = TRUE]
      obj$u <-       obj$u      [, jq, , drop = TRUE]
      obj$k_u <-     obj$k_u    [, jq,   drop = TRUE]
      obj$pi_u <-    obj$pi_u   [, jq, , drop = TRUE]
      if (!is.null(obj$z_u)) {
        obj$z_u <-  obj$z_u[, jq, , drop = TRUE]
      }
      if (length(i) == 1) {
        obj$w <-       matrix(obj$w      , nrow = 1)
        obj$x <-       matrix(obj$x      , nrow = 1)
        obj$theta_x <- matrix(obj$theta_x, nrow = 1)
        obj$theta_v <- matrix(obj$theta_v, nrow = 1)
        obj$theta_e <- matrix(obj$theta_e, nrow = 1)
        obj$var_e <-   matrix(obj$var_e  , nrow = 1)
        obj$x_e <-     matrix(obj$x_e    , nrow = 1)
        obj$x_nz <-    matrix(obj$x_nz   , nrow = 1)
        obj$u <-       matrix(obj$u      , nrow = 1)
        obj$k_u <-     matrix(obj$k_u    , nrow = 1)
        obj$pi_u <-    matrix(obj$pi_u   , nrow = 1)
        if (!is.null(obj$z_u)) {
          obj$z_u <-  matrix(obj$z_u, nrow = 1)
        }
      }
      obj0 <- bdeconvzi_result(unclass(obj)[
        c('w', 'x', 'theta_x', 'theta_v', 'theta_e', 'var_e', 'x_e', 'x_nz',
          'u', 'k_u', 'z_u', 'pi_u', 'p_u', 'mu_u', 'sig1_u', 'sig2_u', 'knots',
          'q', 'p', 'control', 'res_type', 'err_type', 'var_type', 'call')
        ], obj$res_type)
      obj0$q <- length(jq)
      obj0$p <- length(jp)
      obj0$lwr <- obj$lwr[jq]
      obj0$upr <- obj$upr[jq]
      obj0$knots <- obj$knots[jq, ]
      obj0$w_var <- obj$w_var[jq, ]
      obj0$names <- obj$names[jq]
      obj0$y <- obj$y[,jq]
      obj <- obj0
    } else if (length(jp) == 1) {
      obj$w <-       obj$w      [, jp + obj$q, , drop = TRUE]
      obj$x <-       obj$x      [, jp + obj$q, , drop = TRUE]
      obj$z_x <-     obj$z_x    [, jp, ,         drop = TRUE]
      obj$pi_x <-    obj$pi_x   [, jp, ,         drop = TRUE]
      obj$mu_x <-    obj$mu_x   [, jp, ,         drop = TRUE]
      obj$sig_x <-   obj$sig_x  [, jp, ,         drop = TRUE]
      obj$theta_v <- obj$theta_v[, jp + obj$q, , drop = TRUE]
      obj$var_e <-   obj$var_e  [, jp + obj$q,   drop = TRUE]
      obj$u <-       obj$u      [, jp + obj$q, , drop = TRUE]
      obj$k_u <-     obj$k_u    [, jp + obj$q,   drop = TRUE]
      obj$pi_u <-    obj$pi_u   [, jp + obj$q, , drop = TRUE]
      if (!is.null(obj$z_u)) {
        obj$z_u <-  obj$z_u[, jp + obj$q, , drop = TRUE]
      }
      if (length(i) == 1) {
        obj$w <-       matrix(obj$w      , nrow = 1)
        obj$x <-       matrix(obj$x      , nrow = 1)
        obj$pi_x <-    matrix(obj$pi_x   , nrow = 1)
        obj$mu_x <-    matrix(obj$mu_x   , nrow = 1)
        obj$sig_x <-   matrix(obj$sig_x  , nrow = 1)
        obj$theta_v <- matrix(obj$theta_v, nrow = 1)
        obj$var_e <-   matrix(obj$var_e  , nrow = 1)
        obj$u <-       matrix(obj$u      , nrow = 1)
        obj$k_u <-     matrix(obj$k_u    , nrow = 1)
        obj$pi_u <-    matrix(obj$pi_u   , nrow = 1)
        if (!is.null(obj$z_u)) {
          obj$z_u <-  matrix(obj$z_u, nrow = 1)
        }
        if (!is.null(obj$z_x)) {
          obj$z_x <-  matrix(obj$z_x, nrow = 1)
        }
      }
      obj0 <- bdeconvnc_result(unclass(obj)[
        c('w', 'x', 'z_x', 'pi_x', 'mu_x', 'sig_x', 'u', 'k_u', 'z_u', 'pi_u', 'p_u',
          'mu_u', 'sig1_u', 'sig2_u', 'theta_v', 'knots', 'var_e', 'q', 'p',
          'control', 'res_type', 'err_type', 'var_type', 'call')
      ], obj$res_type)
      obj0$q <- length(jq)
      obj0$p <- length(jp)
      obj0$lwr <- obj$lwr[jp + obj$q]
      obj0$upr <- obj$upr[jp + obj$q]
      obj0$knots <- obj$knots[jp + obj$q, ]
      obj0$w_var <- obj$w_var[jp + obj$q, ]
      obj0$names <- obj$names[jp + obj$q]
      obj0$y <- obj$y[,jp + obj$q]
      obj <- obj0
    }
  } else {
    obj$w <-       obj$w[, j, , drop = FALSE]
    obj$x <-       obj$x[, j, , drop = FALSE]
    obj$z_x <-     obj$z_x[, jp, , drop = FALSE]
    obj$pi_x <-    obj$pi_x[, jp, , drop = FALSE]
    obj$mu_x <-    obj$mu_x[, jp, , drop = FALSE]
    obj$sig_x <-   obj$sig_x[, jp, , drop = FALSE]
    obj$theta_x <- obj$theta_x[, jq, , drop = FALSE]
    obj$theta_v <- obj$theta_v[, j, , drop = FALSE]
    obj$x_nz <-    obj$x_nz[, j, , drop = FALSE]
    obj$x_e <-     obj$x_e[, jq, , drop = FALSE]
    obj$theta_e <- obj$theta_e[, jq, , drop = FALSE]
    obj$var_e <-   obj$var_e[, j, drop = FALSE]
    obj$u <-       obj$u[, j, , drop = FALSE]
    obj$k_u <-     obj$k_u[, j, drop = FALSE]
    obj$z_u <-     obj$z_u[, j, , drop = FALSE]
    obj$pi_u <-    obj$pi_u[, j, , drop = FALSE]
    obj$p_u <-     obj$p_u[, , drop = FALSE]
    obj$mu_u <-    obj$mu_u[, , drop = FALSE]
    obj$sig1_u <-  obj$sig1_u[, , drop = FALSE]
    obj$sig2_u <-  obj$sig2_u[, , drop = FALSE]
    obj$corr_x <-  obj$corr_x[, j, j, drop = FALSE]
    obj$corr_e <-  obj$corr_e[, j, j, drop = FALSE]
    obj$q <-       length(jq)
    obj$p <-       length(jp)
    obj$lwr <-     obj$lwr[j]
    obj$upr <-     obj$upr[j]
    obj$knots <-   obj$knots[j, ]
    obj$w_var <-   obj$w_var[j, ]
    obj$names <-   obj$names[j]
    obj$y <- obj$y[,j]
  }
  return(obj)
}
