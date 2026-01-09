#' Error Distribution Evaluation
#'
#' @param obj A `bdeconv_result` object.
#' @param vals A numeric vector (or matrix) on which to evaluate the error
#'  density.
#'
#' @return A numeric matrix or array of cumulative distribution values.
#' @keywords internal
dist_u <- function(obj, vals) {
  UseMethod('dist_u')
}


#' @keywords internal
dist_u.uv_res <- function(obj, vals) {
  if (obj$err_type == 'flex') {
    dist_u <- matrix(0, nrow = length(obj), ncol = length(vals))
    for (i in 1:length(obj)) {
      k_u <- obj[[i]]$k_u
      sd_e <- sqrt(var_e_fn(obj[[i]]$pi_u[1:k_u], obj[[i]]$p_u[1:k_u],
                            obj[[i]]$mu_u[1:k_u], obj[[i]]$sig1_u[1:k_u],
                            obj[[i]]$sig2_u[1:k_u]))
      for (k in 1:k_u) {
        dist_u[i, ] <- dist_u[i, ] + obj[[i]]$pi_u[k] *
          prtcnorm(vals, obj[[i]]$p_u[k], obj[[i]]$mu_u[k] / sd_e,
                   obj[[i]]$sig1_u[k] / sd_e, obj[[i]]$sig2_u[k] / sd_e)
      }
    }
  } else if (obj$err_type == 'norm') {
    dist_u <- matrix(stats::pnorm(vals), nrow = 1, ncol = length(vals))
  } else if (obj$err_type == 'lapl') {
    dist_u <- matrix(ifelse(vals >= 0, exp(vals), 1 - exp(-vals)) / 2,
                     nrow = 1, ncol = length(vals))
  }
  return(dist_u)
}

#' Estimated cumulative distribution function for the true variable(s)
#'
#' @param obj A `bdeconv_result` object.
#' @param vals A numeric vector (for univariate) or matrix (for multivariate) on which to evaluate the CDF.
#' @param expand Used only for `mv_res` objects. A logical value indicating
#'  whether or not to expand the matrix of values.
#'
#' @return A numeric matrix or array of cumulative distribution values for which each row represents the estimated CDF corresponding to each MCMC iteration.
#' @export
dist_x <- function(obj, vals, expand = FALSE) {
  UseMethod('dist_x')
}

#' @export
dist_x.zi_res <- function(obj, vals, expand = FALSE) {
  dist_x <- matrix(0, nrow = length(obj), ncol = length(vals))
  for (i in 1:length(obj)) {
    dist_x[i, ] <- pbspline(vals, obj[[i]]$theta_x, obj$lwr, obj$upr)
  }
  return(dist_x)
}

#' @export
dist_x.nc_res <- function(obj, vals, expand = FALSE) {
  dist_x <- matrix(0, nrow = length(obj), ncol = length(vals))
  for (i in 1:length(obj)) {
    k_x <- length(obj[[i]]$pi_x)
    for (k in 1:k_x) {
      dist_x[i, ] <- dist_x[i, ] + obj[[i]]$pi_x[k] *
        ptnorm(vals, obj[[i]]$mu_x[k], obj[[i]]$sig_x[k], obj$lwr, obj$upr)
    }
  }
  return(dist_x)
}

#' @export
dist_x.mv_res <- function(obj, vals, expand = FALSE) {
  n <- nrow(obj)
  d <- ncol(obj)
  r <- nrow(vals)
  if (expand) {
    idxs <- list()
    for (dd in 1:d) {
      idxs[[dd]] <- 1:r
    }
    grd <- expand.grid(idxs)
  } else {
    grd <- matrix(1:r, nrow = r, ncol = d)
  }
  dens_mvt <- array(0, dim = c(n, nrow(grd)))
  y <- lapply(1:d, function(dd) dist_x(obj[, dd], vals[, dd])) |>
    simplify2array() |>
    pad_probability() |>
    stats::qnorm()
  corr <- obj$corr_x
  for (nn in 1:n) {
    yy <- matrix(0, nrow = nrow(grd), ncol = d)
    for (dd in 1:d) {
      yy[, dd] = y[nn, grd[, dd], dd]
    }
    dens_mvt[nn, ] <- yy |>
      apply(1, function(yyy)
        mvtnorm::pmvnorm(upper = yyy, mean = rep(0,d), corr = corr[nn, , ]))
  }
  return(dens_mvt)
}
