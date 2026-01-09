#' Visualization of `bdeconv_result` object
#' @param x A `bdeconv_result` object.
#' @param ... Unused.
#' @param use_ggplot A logical value indicating whether the plots will be generated using `ggplot2` or using base plot functions.
#' @return A collection of plots illustrating various densities and related functions.
#' @details
#' By default it uses `ggplot2` package if it is available in the system library, otherwise it falls back to the base R plot function.
#'
#' @exportS3Method
plot.bdeconv_result <- function(x, ..., use_ggplot = TRUE) {
  if (requireNamespace('ggplot2', quietly = TRUE) & use_ggplot) {
    p <- autoplot.bdeconv_result(x, ...)
    return(p)
  } else {
    x <- x
    if (length(x) > 100) {
      x <- x[sample(1:length(x), 100)]
    }
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    NextMethod()
  }
}

#' @exportS3Method
plot.uv_res <- function(x, ...) {
  obj <- x
  graphics::par(cex = 0.7, mai = c(0.4, 0.4, 0.1, 0.1), mgp = c(.5, 0.3, 0)) # margin changed to .2 from (0.75, 0.75, 0.1, 0.1), mgp introduced
  len <- length(obj)
  alph <- min(max(5 / len, 0.01), 1)
  graphics::par(fig = c(0.625, 1, 0.45, 0.95))
  value <- seq(-4, 4, length = 500)
  if (obj$err_type == 'flex') {
    dens_u <- dens_u(obj, value)
    if (len > 1){
      plot(value, apply(dens_u, 2, mean), xlab = NA, ylab = NA, xlim = c(-4, 4), type = 'n')
      graphics::polygon(c(value,rev(value)),
                        c(apply(dens_u, 2, stats::quantile, probs = c(.025)),
                          rev(apply(dens_u, 2, stats::quantile, probs = c(.975)))),
                        col = "darkgrey", border = NA)
      graphics::lines(value, stats::dnorm(value), lty = 'dashed')
      graphics::lines(value, apply(dens_u, 2, mean))
    }else{
      plot(value, dens_u, xlab = NA, ylab = NA, type = 'l')
    }
  } else if (obj$err_type == 'norm') {
    plot(value, stats::dnorm(value), xlab = NA, ylab = NA, type = 'l')
  } else if (obj$err_type == 'lapl') {
    plot(value, stats::dexp(abs(value)) / 2, xlab = NA, ylab = NA, type = 'l')
  } else {
    print(obj$err_type)
  }
  graphics::mtext('Scaled Error Value', side = 1, cex = 0.7, line = 1.5)
  graphics::par(fig = c(0.63, 1, 0, 0.45), new = TRUE)
  value <- seq(obj$lwr, obj$upr, length = 100)
  basis <- bspline_basis(value, obj$lwr, obj$upr, length(obj$knots))
  if (len > 1) {
    plot(apply(obj$x, 2, mean), obj$w_var, pch = 3,
         col = grDevices::rgb(0.5, 0.5, 0.5, 0.5), xlab = NA, ylab = NA)
  } else {
    plot(obj$x, obj$w_var, pch = 3,
         col = grDevices::rgb(0.5, 0.5, 0.5, 0.5), xlab = NA, ylab = NA)
  }
  theta_v <- matrix(0, nrow = length(obj), ncol = length(value))
  for (i in 1:len) {
    theta_v[i, ] <- basis %*% t(exp(obj[[i]]$theta_v)) * obj[[i]]$var_e
  }
  graphics::polygon(c(value,rev(value)),
                    c(apply(theta_v, 2, stats::quantile, probs = c(.025)),
                      rev(apply(theta_v, 2, stats::quantile, probs = c(.975)))),
                    col = "darkgrey", border = NA)
  graphics::lines(value, apply(theta_v, 2, mean), col = 'black')
  graphics::mtext('Value', side = 1, cex = 0.7, line = 1.5)
  graphics::mtext('Variance', side = 2, cex = 0.7, line = 1.5)
  NextMethod()
}

#' @exportS3Method
plot.zi_res <- function(x, ...) {
  obj <- x
  len <- length(obj)
  alph <- min(max(5 / len, 0.01), 1)

  graphics::par(fig = c(0, 0.675, 0, 0.35), new = TRUE)
  value0 <- seq(obj$lwr, obj$upr, length = 500)
  prob_z <- matrix(0, nrow = length(obj), ncol = length(value0))
  basis <- bspline_basis(value0, obj$lwr, obj$upr, length(obj[[1]]$theta_v) - 1)
  plot(value0, seq(0, 1, length = length(value0)),
       xlim = c(min(obj$x), max(obj$x)),
       xlab = NA, ylab = NA, type = 'n')
  for (i in 1:len) {
    prob_z[i, ] <- stats::pnorm(basis %*% t(obj[[i]]$theta_e))
  }
  graphics::polygon(c(value0,rev(value0)),
                    c(apply(prob_z, 2, stats::quantile, probs = c(.025)),
                      rev(apply(prob_z, 2, stats::quantile, probs = c(.975)))),
                    col = "darkgrey", border = NA)
  graphics::lines(value0, apply(prob_z, 2, mean), col = 'black')
  graphics::abline(v = obj$lwr, lty = 'dashed')
  graphics::abline(v = obj$upr, lty = 'dashed')
  graphics::mtext('Value', side = 1, cex = 0.7, line = 1.5)
  graphics::mtext('Probability of Zero', side = 2, cex = 0.7, line = 1.5)

  graphics::par(fig = c(0, 0.675, 0.275, 0.95), new = TRUE)
  value <- c(obj$lwr, value0, obj$upr)
  dens_x <- matrix(0, nrow = length(obj), ncol = length(value))
  for (i in 1:len) {
    dens_x[i, ] <- c(0, dbspline(value0, obj[[i]]$theta_x, obj$lwr, obj$upr), 0)
  }
  plot(value, apply(dens_x, 2, mean), xlab= NA, ylab = NA, xlim = c(min(obj$x), max(obj$x)),
                    type = 'n')
  graphics::mtext('Density', side = 2, cex = 0.7, line = 1.5)
  graphics::polygon(c(value,rev(value)),
                    c(apply(dens_x, 2, stats::quantile, probs = c(.025)),
                      rev(apply(dens_x, 2, stats::quantile, probs = c(.975)))),
                    col = "darkgrey", border = NA)
  graphics::lines(value, apply(dens_x, 2, mean))

  graphics::par(cex = 1, mai = c(1.02, 0.82, 0.82, 0.42), fig = c(0, 1, 0, 1))
  graphics::title(main = paste0('Deconvolution of `', obj$names, '`'),
                  line = 3.2)
}

#' @exportS3Method
plot.nc_res <- function(x, ...) {
  obj <- x
  graphics::par(fig = c(0, 0.675, 0, 0.95), new = TRUE)
  len <- length(obj)
  alph <- min(max(5 / len, 0.01), 1)
  value0 <- seq(obj$lwr, obj$upr, length = 500)
  value <- c(obj$lwr, value0, obj$upr)
  dens_x <- matrix(0, nrow = length(obj), ncol = 500)
  for (i in 1:len) {
    k_x <- length(obj[[i]]$pi_x)
    for (k in 1:k_x) {
      dens_x[i, ] <- dens_x[i, ] + obj[[i]]$pi_x[k] *
        dtnorm(value0, obj[[i]]$mu_x[k], obj[[i]]$sig_x[k],
               obj$lwr, obj$upr)
    }
  }
  dens_x <- cbind(0, dens_x, 0)
  plot(value, apply(dens_x, 2, mean),xlab = NA, ylab = NA, xlim = c(min(obj$x), max(obj$x)),
                    type = 'n')
  graphics::polygon(c(value,rev(value)),
                    c(apply(dens_x, 2, stats::quantile, probs = c(.025)),
                      rev(apply(dens_x, 2, stats::quantile, probs = c(.975)))),
                    col = "darkgrey", border = NA)
  graphics::lines(value, apply(dens_x, 2, mean))
  graphics::mtext('Value', side = 1, cex = 0.7, line = 1.5)
  graphics::mtext('Density', side = 2, cex = 0.7, line = 1.5)

  graphics::par(cex = 1, mai = c(1.02, 0.82, 0.82, 0.42), fig = c(0, 1, 0, 1))
  graphics::title(main = paste0('Deconvolution of `', obj$names, '`'),
                  line = 3.2)
}

#' @exportS3Method
plot.mv_res <- function(x, ...) {
  obj <- x
  if (length(obj) > 10) {
    obj <- obj[sample(1:length(obj), 10), ]
  }
  len <- length(obj)
  alph <- min(max(5 / len, 0.01), 1)
  grid_res <- 100
  # Plot univariate results first
  q <- obj$q
  p <- obj$p
  d <- p + q
  var_names <- obj$names
  # Precompute marginal density estimates (can this be passed into plots?)
  dens_x <- array(0, dim = c(len, d, grid_res))
  grid_x <- matrix(0, nrow = d, ncol = grid_res)
  for (nn in 1:len) {
    if (q > 0) {
      for (qq in 1:q) {
        grid_x[qq, ] <- seq(obj$lwr[qq], obj$upr[qq], length = grid_res)
        dens_x[nn, qq, ] <- dbspline(grid_x[qq, ], obj[[nn]]$theta_x[, qq, ],
                                     obj$lwr[qq], obj$upr[qq])
      }
    }
    if (p > 0) {
      for (pp in (q + 1:p)) {
        grid_x[pp, ] <- seq(obj$lwr[pp], obj$upr[pp], length = grid_res)
        k_x <- length(obj[[nn]]$pi_x[, pp - q, ])
        for (kk in 1:k_x) {
          dens_x[nn, pp, ] <- dens_x[nn, pp, ] + obj[[nn]]$pi_x[, pp - q, kk] *
            dtnorm(grid_x[pp, ], obj[[nn]]$mu_x[, pp - q, kk],
                   obj[[nn]]$sig_x[, pp - q, kk], obj$lwr[pp], obj$upr[pp])
        }
      }
    }
  }
  # Compute joint density estimates
  grid_y <- array(0, dim = c(len, d, grid_res))
  dist_x <- array(0, dim = c(len, d, grid_res))
  for (dd in 1:d) {
    for (nn in 1:len) {
      dist_x[nn, dd, ] <- cumsum(dens_x[nn, dd, ]) *
        (grid_x[dd, 2] - grid_x[dd, 1])
      dist_x[nn, dd, dist_x[nn, dd, ] > 1] <- 1
      grid_y[nn, dd, ] <- stats::qnorm(dist_x[nn, dd, ])
    }
  }
  grid_y[grid_y == Inf] <- 3.65
  grid_y[grid_y == -Inf] <- -3.65

  idxs <- utils::combn(d, 2)
  grid_y2 <- array(0, dim = c(2, grid_res ^ 2))
  dens_x2 <- array(0, dim = c(2, grid_res ^ 2))
  dens_x2arr <- array(0, dim = c(len, d, d, grid_res, grid_res))
  for (nn in 1:len) {
    for (kk in 1:ncol(idxs)) {
      dim_idxs <- idxs[, kk]
      d1 <- dim_idxs[1]
      d2 <- dim_idxs[2]

      dens_x2[1, ] <- rep(dens_x[nn, d1, ], times = grid_res)
      dens_x2[2, ] <- rep(dens_x[nn, d2, ], each = grid_res)

      grid_y2 <- rbind(rep(grid_y[nn, d1, ], time = grid_res),
                       rep(grid_y[nn, d2, ], each = grid_res))

      dens_x2arrvec <- (dmvnorm(t(grid_y2), c(0, 0),
                                obj[[nn]]$corr_x[, c(d1 ,d2), c(d1, d2)], FALSE)
                        / dmvnorm(t(grid_y2), c(0, 0), diag(1, 2), FALSE)) *
        apply(dens_x2, 2, prod)
      dens_x2arrvec[is.na(dens_x2arrvec)] <- 0
      dens_x2arr[nn, d1, d2, , ] <- matrix(dens_x2arrvec, nrow = grid_res,
                                           byrow = TRUE)
    }
  }
  # Plot joint results
  laymat <- diag(1:d)
  laymat[lower.tri(laymat)] <- (d + 1):(d * (d + 1) / 2)
  laymat <- t(laymat)
  laymat_tmp <- diag(rep(0, d))
  laymat_tmp[upper.tri(laymat_tmp)] <- ((d * (d + 1) / 2) + 1):d ^ 2
  laymat <- laymat + t(laymat_tmp)
  graphics::layout(laymat)

  graphics::par(mar = rep(1.85, 4))
  for (dd in 1:d)  {
    plot(grid_x[dd, ], apply(dens_x[, dd, , drop = FALSE], 3, mean),
                      ylim = c(0, max(dens_x[, dd, ])),
                      type = 'n', xlab = NA, ylab = NA)
    graphics::polygon(c(grid_x[dd, ],rev(grid_x[dd, ])),
                      c(apply(dens_x[, dd, , drop = FALSE], 3, stats::quantile, probs = c(.025)),
                        rev(apply(dens_x[, dd, , drop = FALSE], 3, stats::quantile, probs = c(.975)))),
                      col = "darkgrey", border = NA)
    graphics::lines(grid_x[dd, ], apply(dens_x[, dd, , drop = FALSE], 3, mean))
    if (dd == 1){
      graphics::mtext(var_names[dd], side = 3, line = .5, cex = 0.8)
    }
    if (dd == d){
      graphics::mtext(var_names[d], side = 4, line = .5, las = 3, cex = 0.8)
    }
  }
  grid_lim <- 1:100 * grid_res / 100
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      for (nn in 1:len) {
        graphics::contour(grid_x[j, grid_lim], grid_x[i, grid_lim],
                          dens_x2arr[nn, i, j, grid_lim, grid_lim],
                          add = (nn != 1), nlevels = 15, drawlabels = FALSE,
                          col = grDevices::rgb(0, 0, 0, alph / 5))
      }
      graphics::contour(grid_x[j, grid_lim], grid_x[i, grid_lim],
                        apply(
                          dens_x2arr[, i, j, grid_lim, grid_lim, drop = FALSE],
                          4:5,
                          mean
                        ),
                        add = TRUE, nlevels = 15, drawlabels = FALSE,
                        col = 'black')
      if (i==1){
        graphics::mtext(var_names[j], side = 3, line = .5, cex = 0.8)
      }
      if (j==d){
        graphics::mtext(var_names[i], side = 4, line = .5, las = 3, cex = 0.8)
      }
    }
  }
  for (j in 1:(d - 1)) {
    for (i in (j + 1):d) {
      plot(apply(obj$x[, i, , drop = FALSE], 3, mean),
           apply(obj$x[, j, , drop = FALSE], 3, mean),
           pch = 3, col = grDevices::rgb(0.5, 0.5, 0.5, 0.25),
           xlab = NA, ylab = NA)
    }
  }
}

