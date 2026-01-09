#' Visualization of `bdeconv_result` object using `ggplot2`
#' @param object A `bdeconv_result` object.
#' @return A collection of plots illustrating various densities and related functions made using `ggplot2`.
#' @export
autoplot.bdeconv_result <- function(object) {
  if (length(object) > 100) {
    object <- object[sample(1:length(object), 100)]
  }
  if (inherits(object, 'uv_res')) {
    p1 <- .autoplot_errdens(object)
    p2 <- .autoplot_errfunc(object)
    p3 <- .autoplot_valdens(object)
    if (inherits(object, 'zi_res')) {
      p4 <- .autoplot_zerprob(object)
      gridExtra::grid.arrange(gridExtra::arrangeGrob(p3, p4,
                                                     ncol = 1,
                                                     heights = c(2, 1)),
                              gridExtra::arrangeGrob(p1, p2,
                                                     ncol = 1,
                                                     heights = c(1, 1)),
                              nrow = 1,
                              top = paste0('Density Deconvolution for `',
                                           object$name, '`'))
    } else if (inherits(object, 'nc_res')) {
      gridExtra::grid.arrange(p3,
                              gridExtra::arrangeGrob(p1, p2,
                                                     ncol = 1,
                                                     heights = c(1, 1)),
                              nrow = 1,
                              top = paste0('Density Deconvolution for `',
                                           object$name, '`'))
    }
  } else if (inherits(object, 'mv_res')) {
    d <- dim(object)[2]
    ps <- list()
    for (i in 1:d) {
      for (j in 1:d) {
        ps <- append(ps, list(.autoplot_mvgrid(object, i, j)))
      }
    }
    gridExtra::grid.arrange(grobs = ps, nrow = d,
                            top = paste('Density Deconvolution for',
                                        paste(
                                          paste0('`', object$names, '`'),
                                          collapse = ', '),
                                        collapse = ' '))
  }
}

.autoplot_errdens <- function(obj) {
  vals <- seq(-4, 4, length = 500)
  dens_u <- dens_u(obj, vals)
  if ((obj$err_type == 'flex') && (length(obj) > 1)) {

    df <- data.frame(Iteration = rep(1, length(vals)),
                     Value     = vals,
                     Density   = apply(dens_u,2,mean))
    alpha <- 1
    color <- 'black'
  } else {
    df <- data.frame(Iteration = rep(1, length(vals)),
                     Value     = vals,
                     Density   = as.numeric(dens_u))
    alpha <- 1
    color <- 'black'
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Value"]], y = .data[["Density"]],
                                               group = .data[["Iteration"]]))
  if (obj$err_type == 'flex') {
    if (length(obj) > 1) {
      p <- p + ggplot2::geom_ribbon(
        data = data.frame(value = vals,
                          lower=apply(dens_u, 2, stats::quantile, probs = 0.025),
                          upper=apply(dens_u, 2, stats::quantile, probs = 0.975)),
        ggplot2::aes(x = .data[["value"]], ymin = .data[["lower"]], ymax = .data[["upper"]]),
        fill = 'grey', alpha = 1, inherit.aes = FALSE)
    }
    p <- p + ggplot2::geom_line(data = data.frame(Iteration = rep(0,length(vals)),
                                                  Value = vals,
                                                  Density = stats::dnorm(vals)),
                                linetype = 'dashed')
  }
  p <- p + ggplot2::geom_line(color = color, alpha = alpha) +
    ggplot2::xlab("Scaled Error Value")+
    ggplot2::theme_bw()
  return(p)
}

.autoplot_errfunc <- function(obj) {
  alpha <- min(max(5 / length(obj), 0.01), 1)
  color <- ifelse(length(obj) == 1, 'black', 'grey')
  x <- apply(obj$x, 2, mean)
  vals <- seq(min(x), max(x), length = 500)
  var_u <- varfun(obj, vals)
  df <- data.frame(Iteration = rep(1:length(obj), length(vals)),
                   Value     = rep(vals, each = length(obj)),
                   Variance  = as.numeric(var_u))
  df <- df[df$Variance <= max(obj$w_var, na.rm = TRUE), ]
  df_x <- data.frame(Iteration = 0,
                     Value    = x,
                     Variance = obj$w_var)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Value"]], y = .data[["Variance"]],
                                               group = .data[["Iteration"]]))+
    ggplot2::theme_bw()
  if (length(obj) > 1) {
    df_q <- data.frame(value = vals,
                       lower=apply(var_u, 2, stats::quantile, probs = 0.025),
                       upper=apply(var_u, 2, stats::quantile, probs = 0.975))
    p <- p + ggplot2::geom_ribbon(
      data = stats::na.omit(df_q),
      ggplot2::aes(x = .data[["value"]], ymin = .data[["lower"]], ymax = .data[["upper"]]),
      fill = 'grey', alpha = 1, inherit.aes = FALSE)

    df_m <- data.frame(Iteration = 0,
                       Value     = vals,
                       Variance  = apply(var_u, 2, mean))
    p <- p + ggplot2::geom_line(data = stats::na.omit(df_m))
  } else {
    p <- p + ggplot2::geom_line(color = color, alpha = alpha)
  }
  return(p)
}

.autoplot_valdens <- function(obj) {
  vals <- seq(obj$lwr, obj$upr, length = 1500)
  dens_x <- dens_x(obj, vals)
  vals <- c(obj$lwr, vals, obj$upr)
  dens_x <- cbind(0, dens_x, 0)
  if (length(obj) > 1) {
    df <- data.frame(Iteration = rep(1, length(vals)),
                     Value     = vals,
                     Density   = apply(dens_x, 2, mean))
    alpha <- 1
    color <- 'black'
  } else {
    df <- data.frame(Iteration = rep(1, length(vals)),
                     Value     = vals,
                     Density   = as.numeric(dens_x))
    alpha <- 1
    color <- 'black'
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Value"]], y = .data[["Density"]],
                                               group = .data[["Iteration"]]))
  if (length(obj) > 1) {

    p <- p + ggplot2::geom_ribbon(
      data = data.frame(value = vals,
                        lower=apply(dens_x, 2, stats::quantile, probs = 0.025),
                        upper=apply(dens_x, 2, stats::quantile, probs = 0.975)),
      ggplot2::aes(x = .data[["value"]], ymin = .data[["lower"]], ymax = .data[["upper"]]),
      fill = 'grey', alpha = 1, inherit.aes = FALSE)
  }
  p <- p + ggplot2::geom_line(color = color, alpha = alpha) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = obj$name, y = "density")
  return(p)
}

.autoplot_zerprob <- function(obj) {
  vals <- seq(obj$lwr, obj$upr, length = 500)
  prob_nz <- matrix(0, nrow = length(obj), ncol = length(vals))
  basis <- bspline_basis(vals, obj$lwr, obj$upr, length(obj[[1]]$theta_e) - 1)
  for (i in 1:length(obj)) {
    prob_nz[i, ] <- stats::pnorm(basis %*% t(obj[[i]]$theta_e))
  }
  if (length(obj) > 1) {

    df <- data.frame(Iteration    = rep(1, length(vals)),
                     Value        = vals,
                     Probability  = apply(prob_nz, 2, mean))
    alpha <- 1
    color <- 'black'
  } else {
    df <- data.frame(Iteration = rep(1, length(vals)),
                     Value     = vals,
                     Probability   = as.numeric(prob_nz))
    alpha <- 1
    color <- 'black'
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Value"]], y = .data[["Probability"]],
                                               group = .data[["Iteration"]]))
  if (length(obj) > 1) {
    p <- p + ggplot2::geom_ribbon(
      data = data.frame(value = vals,
                        lower=apply(prob_nz, 2, stats::quantile, probs = 0.025),
                        upper=apply(prob_nz, 2, stats::quantile, probs = 0.975)),
      ggplot2::aes(x = .data[["value"]], ymin = .data[["lower"]], ymax = .data[["upper"]]),
      fill = 'grey', alpha = 1, inherit.aes = FALSE)
  }
  p <- p + ggplot2::geom_line(color = color, alpha = alpha) +
    ggplot2::ylab("Probability of zero")+
    ggplot2::theme_bw()
  return(p)
}

.autoplot_mvgrid <- function(obj, i, j, res = 100, fun = c('dens', 'dist')) {
  if (!inherits(obj, 'mv_res')) {
    stop('`obj` must be an `mv_res` object.')
  }
  fun <- match.arg(fun)
  if (i == j) {
    return(.autoplot_valdens(obj[, i]))
  } else if (i < j) {
    return(.autoplot_contour(obj, i, j, res = res, fun = fun))
  } else if (i > j) {
    return(.autoplot_scatter(obj, i, j))
  }
}

.autoplot_contour <- function(obj, i, j, res = 100, fun = c('dens', 'dist')) {
  if (!inherits(obj, 'mv_res')) {
    stop('`obj` must be an `mv_res` object.')
  }
  if ((i >= j)) {
    stop('`i` must be less than `j`.')
  }
  if (any(c(i, j) <= 0)) {
    stop('`i` and `j` must be positive.')
  }
  fun <- match.arg(fun)

  #obj <- mean(obj)
  vals <- sapply(c(i, j), function(dd)
    seq(obj$lwr[dd], obj$upr[dd], length = res + 2)
    )
  # Remove evaluation at lwr and upr to avoid numerical errors at boundaries
  vals <- vals[2:(res - 1), ]
  if (fun == 'dens') {
    dens <- colMeans(dens_x(obj[, c(i, j)], vals, expand = TRUE), na.rm = TRUE)
  } else if (fun == 'dist') {
    dens <- colMeans(dist_x(obj[, c(i, j)], vals, expand = TRUE), na.rm = TRUE)
  }

  df <- expand.grid(vals[, 1], vals[, 2])
  df$Density <- dens

  my.cols <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  contour.cols <- my.cols[6:11]
  contour_breaks <- pretty(range(df$Density, na.rm = TRUE), n = 45)
  names(df)[1:2] <- obj$names[c(i, j)]
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[obj$names[j]]],
                                          y = .data[[obj$names[i]]],
                                          z = .data[["Density"]])) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data[["Density"]]))+
    ggplot2::geom_contour(breaks = contour_breaks,
                          ggplot2::aes(color = as.factor(ggplot2::after_stat(level))),
                          linewidth = 0.6,
                          show.legend = FALSE) +
    ggplot2::scale_color_manual(values = rep(contour.cols, length.out = length(contour_breaks))) +
    ggplot2::scale_fill_gradientn(colours = grDevices::terrain.colors(100)) +
    ggplot2::guides(fill = 'none') +
    ggplot2::theme_bw() +
    ggplot2::labs(x = obj$names[j], y = obj$names[i])
}

.autoplot_scatter <- function(obj, i, j) {
  if (!inherits(obj, 'mv_res')) {
    stop('`obj` must be an `mv_res` object.')
  }
  if ((i <= j)) {
    stop('`i` must be greater than `j`.')
  }
  if (any(c(i, j) <= 0)) {
    stop('`i` and `j` must be positive.')
  }

  obj <- mean(obj)
  df <- data.frame(xi = obj$x[, i, ],
                   xj = obj$x[, j, ])
  names(df)[1:2] <- obj$names[c(i, j)]
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[obj$names[i]]], y = .data[[obj$names[j]]])) +
    ggplot2::geom_point(shape = 3, color = 'grey50', alpha = 0.25) +
    ggplot2::theme_bw()
}

#' Register `autoplot` methods to `ggplot2`
#'
#' @return None
#' @keywords internal
register_autoplot_s3_methods = function() {
  if (requireNamespace('ggplot2', quietly = TRUE)) {
    register_s3_method('ggplot2', 'autoplot', 'bdeconv_result')
  }
}
