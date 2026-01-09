#' Error Density Evaluation
#'
#' @param obj A `bdeconv_result` object.
#' @param vals A numeric vector (or matrix) on which to evaluate the error
#'  density.
#'
#' @return A numeric matrix or array of density values.
#' @keywords internal
dens_u <- function(obj, vals) {
  UseMethod('dens_u')
}

#' @keywords internal
dens_u.uv_res <- function(obj, vals) {
  if (obj$err_type == 'flex') {
    dens_u <- matrix(0, nrow = length(obj), ncol = length(vals))
    for (i in 1:length(obj)) {
      k_u <- obj[[i]]$k_u
      sd_e <- sqrt(var_e_fn(obj[[i]]$pi_u[1:k_u], obj[[i]]$p_u[1:k_u],
                            obj[[i]]$mu_u[1:k_u], obj[[i]]$sig1_u[1:k_u],
                            obj[[i]]$sig2_u[1:k_u]))
      for (k in 1:k_u) {
        dens_u[i, ] <- dens_u[i, ] + obj[[i]]$pi_u[k] *
          drtcnorm(vals, obj[[i]]$p_u[k], obj[[i]]$mu_u[k] / sd_e,
                   obj[[i]]$sig1_u[k] / sd_e, obj[[i]]$sig2_u[k] / sd_e)
      }
    }
  } else if (obj$err_type == 'norm') {
    dens_u <- matrix(stats::dnorm(vals), nrow = 1, ncol = length(vals))
  } else if (obj$err_type == 'lapl') {
    dens_u <- matrix(stats::dexp(abs(vals)) / 2, nrow = 1, ncol = length(vals))
  }
  return(dens_u)
}

#' @keywords internal
dens_u.mv_res <- function(obj, vals) {
  return(dens_mv(obj, vals, type = 'u'))
}

#' Estimated density function for the true variable(s)
#'
#' @param obj A `bdeconv_result` object.
#' @param vals A numeric vector (for univariate) or matrix (for multivariate) on which to evaluate the density.
#' @param expand Used only for `mv_res` objects. A logical value indicating
#'  whether or not to expand the matrix of values.
#'
#'
#' @return A numeric matrix or array of density values for which each row represents the estimated density corresponding to each MCMC iteration.
#' @export
dens_x <- function(obj, vals, expand = FALSE) {
  UseMethod('dens_x')
}

#' @export
dens_x.zi_res <- function(obj, vals, expand = FALSE) {
  dens_x <- matrix(0, nrow = length(obj), ncol = length(vals))
  for (i in 1:length(obj)) {
    dens_x[i, ] <- dbspline(vals, obj[[i]]$theta_x, obj$lwr, obj$upr)
  }
  return(dens_x)
}

#' @export
dens_x.nc_res <- function(obj, vals, expand = FALSE) {
  dens_x <- matrix(0, nrow = length(obj), ncol = length(vals))
  for (i in 1:length(obj)) {
    k_x <- length(obj[[i]]$pi_x)
    for (k in 1:k_x) {
      tmp <- obj[[i]]$pi_x[k] *
        dtnorm(vals, obj[[i]]$mu_x[k], obj[[i]]$sig_x[k], obj$lwr, obj$upr)
      tmp <- replace(tmp, is.nan(tmp), 0)
      tmp <- replace(tmp, is.infinite(tmp), 0)
      dens_x[i, ] <- dens_x[i, ] + replace(tmp, is.nan(tmp), 0)
    }
  }
  return(dens_x)
}

#' @export
dens_x.mv_res <- function(obj, vals, expand = FALSE) {
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
  dens_uni <- lapply(1:d, function(dd) dens_x(obj[, dd], vals[, dd])) |>
    simplify2array()
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
    prec <- solve(corr[nn, , ]) - diag(d)
    dens_mvt[nn, ] <- (dmvnorm(yy, rep(0,d), corr[nn, , ], FALSE)
                      / dmvnorm(yy, rep(0,d), diag(1, d), FALSE)) *
      apply(
        sapply(1:d, function(dd)
          dens_uni[nn, grd[, dd], dd]),
        1, prod)
  }
  return(dens_mvt)
}


dens_x_normed <- function(obj, vals, norm_by) {
  obj <- mean(obj)
  n <- dim(obj$x)[[3]]
  d <- ncol(obj)
  r <- nrow(vals)
  dists <- grid_y <- matrix(0, nrow = d, ncol = r)
  for(dd in 1:d) {
    dists[dd, ] <- dist_x(obj[, dd], vals[, dd])
    grid_y[dd, ] <- stats::qnorm(dists[dd, ])
  }
  x_norm <- matrix(1, nrow = d, ncol = n)
  x_grid <- matrix(1, nrow = d, ncol = r)
  for(dd in 1:d) {
    if (dd != norm_by) {
      x_norm[dd, ] <- obj$x[, dd, ] / obj$x[, norm_by, ]
      x_grid[dd,] <- seq(0, stats::quantile(x_norm[dd, ], 0.9999), len = r)
    }
  }
  dens_norm <- matrix(0, nrow = d - 1, ncol = r)
  y_norm <- matrix(0, nrow = d, ncol = r)
  y_norm[norm_by, ] <- grid_y[norm_by, ]
  denss <- t(sapply(1:d, function(dd) dens_x(obj[, dd],vals[, dd])))
  dens_norm <- dens_norm_tmp <- denss
  for(dd in 1:d) {
    if (dd != norm_by) {
      for(rr in 1:r) {
        x.dd.grid.rr <- x_grid[dd, rr] * vals[, norm_by]
        TempMat <- abs(matrix(rep(x.dd.grid.rr, r), r, r) -
                         matrix(rep(vals[, dd], n), r, r, byrow = TRUE))
        close.ind <- apply(TempMat, 1, which.min)

        y_norm[dd, ] <- grid_y[dd, close.ind]
        dens_norm_tmp[dd, ] <- denss[dd, close.ind]
        dens_biv_norm <- (dmvnorm(t(y_norm[c(dd, norm_by), ]), rep(0, 2),
                                  obj$corr_x[, c(dd, norm_by), c(dd, norm_by)],
                                  FALSE) /
                            dmvnorm(t(y_norm[c(dd, norm_by),] ), rep(0,2),
                                    diag(1,2), FALSE)) *
          apply(dens_norm_tmp[c(dd, norm_by), ], 2, prod)
        dens_biv_norm[which(is.nan(dens_biv_norm))] <- 0
        dens_norm[dd, rr] <- sum(vals[, dd] * dens_biv_norm) *
          (vals[2, dd] - vals[1, dd])
      }
    }
  }
  return(dens_norm)
}

dens_mv <- function(obj, vals, type = c('x', 'u'), norm_by = 0) {
  type <- match.arg(type)
  n <- nrow(obj)
  d <- ncol(obj)
  r <- nrow(vals)
  dens_v2arr <- array(0, dim = c(n, d, d, r, r))
  for (nn in 1:n) {
    if (type == 'x') {
      dists <- t(sapply(1:d, function(dd) dist_u(obj[nn, dd], vals[, dd])))
      denss <- t(sapply(1:d, function(dd) dens_u(obj[nn, dd], vals[, dd])))
    } else if (type == 'u') {
      dists <- t(sapply(1:d, function(dd) dist_u(obj[nn, dd], vals[, dd])))
      denss <- t(sapply(1:d, function(dd) dens_u(obj[nn, dd], vals[, dd])))
    }
    grid_y <- apply(dists, 1, stats::qnorm) |> t()
    idxs <- utils::combn(d, 2)
    grid_y2 <- array(0, dim = c(2, r ^ 2))
    dens_v2 <- array(0, dim = c(2, r ^ 2))
    dens_v2arr <- array(0, dim = c(n, d, d, r, r))

    for (kk in 1:ncol(idxs)) {
      dim_idxs <- idxs[, kk]
      d1 <- dim_idxs[1]
      d2 <- dim_idxs[2]

      dens_v2[1, ] <- rep(denss[d1, ], times = r)
      dens_v2[2, ] <- rep(denss[d2, ], each = r)
      grid_y2 <- rbind(rep(grid_y[d1, ], times = r),
                       rep(grid_y[d2, ], each = r))
      if (type == 'x') {
        corr <- obj[[nn]]$corr_x[, c(d1 ,d2), c(d1, d2)]
      } else if (type == 'u') {
        corr <- obj[[nn]]$corr_e[, c(d1 ,d2), c(d1, d2)]
      }

      dens_v2arrvec <- (dmvnorm(t(grid_y2), c(0, 0), corr, FALSE)
                        / dmvnorm(t(grid_y2), c(0, 0), diag(1, 2), FALSE)) *
        apply(dens_v2, 2, prod)
      dens_v2arrvec[is.na(dens_v2arrvec)] <- 0
      dens_v2arr[nn, d1, d2, , ] <- matrix(dens_v2arrvec, nrow = r,
                                           byrow = TRUE)
    }
  }
  return(dens_v2arr)
}

pad_probability <- function(x) {
  x[x <= stats::pnorm(-3.5)] <- stats::pnorm(-3.5)
  x[x >= stats::pnorm(3.5)]  <- stats::pnorm(3.5)
  return(x)
}

pad_quantile <- function(y) {
  y[y == -Inf] <- -3.65
  y[y == Inf] <- 3.65
  return(y)
}
