#' Maximization function for variance function B-spline parameter initialization
#'  in [bdeconv()]
#'
#' @param thetas Spline parameters.
#' @param xs Data values.
#' @param mis Subject indicators.
#' @param knots.t Knot points.
#' @param P.t Penalty matrix.
#' @param s2t Variance.
#' @param us Residuals.
#'
#' @return A score to be maximized.
#' @keywords internal
fr <- function(thetas, xs, mis, knots.t, P.t, s2t, us) {
  vars <- bspline_basis(xs, min(knots.t), max(knots.t), length(knots.t)) %*%
    exp(thetas)
  y = -t(thetas) %*% P.t %*% thetas / (2 * s2t) -
    sum(rep(log(vars), times = mis)) / 2 -
    sum(us ^ 2 / rep(vars, times = mis)) / 2
  return(-y)
}

#' Maximization function for B-spline distribution parameter initialization in
#'  [bdeconv()]
#'
#' @param thetas Spline parameters.
#' @param xs Data values.
#' @param knots.t Knot points.
#' @param P.t Penalty matrix.
#' @param s2t Variance
#' @param density.est Density estimate on grid.
#' @param K.t Spline cardinality.
#'
#' @return A score to be maximized.
#' @keywords internal
fr.dens <- function(thetas, xs, knots.t, P.t, s2t, density.est, K.t) {
  delta = knots.t[2] - knots.t[1]
  density.est2 = bspline_basis(xs, min(knots.t), max(knots.t), K.t) %*%
    bspline_normcoef(thetas, min(knots.t), max(knots.t))
  y = -t(thetas) %*% P.t %*% thetas / (2 * s2t) -
    sum((density.est - density.est2) ^ 2)
  return(-y)
}

#' Maximization function for B-spline zero-inflation probability function
#'  parameter initialization in [bdeconv()]
#'
#' @param log_diff_thetas Spline parameters.
#' @param xs Data values.
#' @param knots.t Knot points.
#' @param P.t Penalty matrix.
#' @param s2t Variance
#' @param probs.emp Empirical consumption probability.
#'
#' @return A score to be maximized.
#' @keywords internal
fr.prob.nonzero <- function(log_diff_thetas, xs, knots.t, P.t, s2t, probs.emp) {
  thetas = cumsum(exp(log_diff_thetas))
  probs.est = as.vector(stats::pnorm(bspline_basis(xs, min(knots.t),
                                                   max(knots.t),
                                                   length(knots.t)) %*% thetas))
  y = -t(thetas) %*% P.t %*% thetas / (2 * s2t) -
    sum((probs.emp - probs.est) ^ 2)
  return(-y)
}



fr.dtnorm <- function(params.xs, xs, lwr, upr){
  mu_ <- params.xs[1]
  sig_ <- exp(params.xs[2])
  loglik <- stats::dnorm(xs,mu_,sig_, log = TRUE) - log((stats::pnorm(upr,mu_,sig_) - stats::pnorm(lwr,mu_,sig_)))
  return(-sum(loglik))
}
