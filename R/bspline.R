#' The Quadratic B-Spline Distribution
#'
#' Density function for the normalized quadratic B-Spline distribution.
#'
#' @aliases bspline
#' @aliases dbspline
#'
#' @name B-Spline
#'
#' @param x A vector of quantiles.
#' @param theta A vector of unnormalized spline parameters.
#' @param lwr A scalar lower bound.
#' @param upr A scalar upper bound.
#' @param log A logical value indication whether a probability `p` should be
#'  given as `log(p)`
#'
#' @return `dbspline` gives the density. The length of the result is determined
#'  by the length of `x` for `dbspline`.
#'
#' @rdname bspline
#' @keywords internal
dbspline <- function(x, theta, lwr, upr, log = FALSE) {
  stopifnot(lwr < upr)
  stopifnot(length(theta) > 5)

  K <- length(theta) - 1
  y <- numeric(length(x))
  idx <- which(x >= lwr & x <= upr)
  y[-idx] <- 0
  basis <- bspline_basis(x[idx], lwr, upr, K)
  y[idx] <- basis %*% bspline_normcoef(theta, lwr, upr)
  if (log) {
    return(log(y))
  } else {
    return(y)
  }
}

#' @rdname bspline
#' @keywords internal
pbspline <- function(x, theta, lwr, upr, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(lwr < upr)
  stopifnot(length(theta) > 5)

  K <- length(theta) - 1
  y <- numeric(length(x))
  idx <- which(x >= lwr & x <= upr)
  y[x < lwr] <- 0
  y[x > upr] <- 1
  basis <- bspline_cumbasis(x[idx], lwr, upr, K)
  y[idx] <- basis %*% bspline_normcoef(theta, lwr, upr)
  if (log.p) {
    return(log(y))
  } else {
    return(y)
  }
}

#' Generate a B-spline Basis Matrix
#'
#' @param x A numeric vector of values on which to evaluate the basis.
#' @param lwr A lower bound.
#' @param upr An upper bound.
#' @param K A number of knot points.
#'
#' @return A B-spline basis matrix.
#' @seealso [bspline]
#' @keywords internal
bspline_basis <- function(x, lwr, upr, K) {
  n <- length(x)
  knots <- seq(lwr, upr, length = K)
  B <- matrix(0, n, K + 1)
  for (jj in 1:(K - 1)) {
    act.inds <- (1:n)[(x >= knots[jj]) & (x <= knots[jj + 1])]
    act.x <- x[act.inds]
    resc.x <- (act.x - knots[jj]) / (knots[jj + 1] - knots[jj])

    B[act.inds, jj] <- (1 / 2) * (1 - resc.x) ^ 2
    B[act.inds, jj + 1] <- -(resc.x ^ 2) + resc.x + 1 / 2
    B[act.inds, jj + 2] <- (resc.x ^ 2) / 2
  }
  return(B)
}

#' Generate a B-spline Cumulative Basis Matrix
#'
#' @param x A numeric vector of values on which to evaluate the basis.
#' @param lwr A lower bound.
#' @param upr An upper bound.
#' @param K A number of knot points.
#'
#' @return A B-spline cumulative basis matrix.
#' @seealso [bspline]
#' @keywords internal
bspline_cumbasis <- function(x, lwr, upr, K) {
  n <- length(x)
  knots <- seq(lwr, upr, length = K)
  delta <- knots[2] - knots[1]
  B <- matrix(0, n, K + 1)
  for (jj in 1:(K-1)) {
    act.inds <- (1:n)[(x >= knots[jj])]
    act.x <- x[act.inds]
    resc.x <- (act.x - knots[jj]) / (knots[jj + 1] - knots[jj])
    resc.x[resc.x > 1] <- 1

    B[act.inds, jj] <- 5 * delta / 6 +
      (delta / 2) * (resc.x - resc.x ^ 2 + resc.x ^ 3 / 3)
    B[act.inds, jj + 1] <- delta / 6 +
      delta * (-resc.x ^ 3 / 3 + resc.x ^ 2 / 2 + resc.x / 2)
    B[act.inds, jj + 2] <- delta * resc.x ^ 3 / 6
  }
  act.inds <- (1:n)[(x >= knots[1])]
  B[act.inds, 1] <- B[act.inds, 1] - 5 * delta / 6
  B[act.inds, 2] <- B[act.inds, 2] - delta / 6
  return(B)
}

#' Normalize B-Spline Parameters
#'
#' @param theta A B-spline parameter vector.
#' @param lwr A lower bound.
#' @param upr An upper bound.
#'
#' @return Normalize B-spline parameters.
#' @seealso [bspline]
#' @keywords internal
bspline_normcoef <- function(theta, lwr, upr) {
  K <- length(theta) - 1
  delta <- (upr - lwr) / (K - 1)
  coeffs <- exp(theta) /
    ((delta / 6) * (exp(theta[1]) +
     5 * exp(theta[2]) +
     6 * sum(exp(theta[3:(K - 1)])) +
     5 * exp(theta[K]) +
     exp(theta[K + 1])))
  return(as.numeric(coeffs))
}
