#' The Restricted Two-Component Normal Distribution
#'
#' Density and distribution functions for the mean-restricted
#'  two-component normal distribution.
#'
#' @aliases rtcnorm
#' @aliases drtcnorm
#' @aliases prtcnorm
#'
#' @name RTCNormal
#'
#' @param x A vector of quantiles.
#' @param p A scalar mixture weight.
#' @param mu A scalar separation parameter.
#' @param sig1 A scalar standard deviation for the first component.
#' @param sig2 A scalar standard deviation for the second component.
#' @param lower.tail A logical value indicated whether or not to take the
#'  cumulative distribution function with respect to the lower tail
#'   (alternatively upper).
#' @param log A logical value indication whether a probability `p`.
#'
#' @return `drtcnorm` gives the density, `prtcnorm` gives the distribution.
#'  The length of the result is determined by the length of `x` for `drtcnorm`.

#' @rdname rtcnorm
#' @keywords internal
drtcnorm <- function(x, p, mu, sig1, sig2, log = FALSE) {
  if (!all(sapply(list(p, mu, sig1, sig1), length) == 1)) {
    stop('`mu`, `sig1`, and `sig2` must have length 1')
  }
  c0 <- p ^ 2 + (1 - p) ^ 2
  c1 <- (1 - p) / c0
  c2 <- -p / c0
  mu1 <- c1 * mu
  mu2 <- c2 * mu
  d <- p * stats::dnorm(x, mu1, sig1) +
    (1 - p) * stats::dnorm(x, mu2, sig2)
  if (log) {
    return(log(d))
  } else {
    return(d)
  }
}

#' @rdname rtcnorm
#' @keywords internal
prtcnorm <- function(q, p, mu, sig1, sig2, lower.tail = TRUE, log.p = FALSE) {
  if (!all(sapply(list(p, mu, sig1, sig1), length) == 1)) {
    stop('`mu`, `sig1`, and `sig2` must have length 1')
  }
  c0 <- p ^ 2 + (1 - p) ^ 2
  c1 <- (1 - p) / c0
  c2 <- -p / c0
  mu1 <- c1 * mu
  mu2 <- c2 * mu
  y <- p * stats::pnorm(q, mu1, sig1, lower.tail = lower.tail) +
    (1 - p) * stats::pnorm(q, mu2, sig2, lower.tail = lower.tail)
  if (log.p) {
    return(log(y))
  } else {
    return(y)
  }
}

#' Variance Scaling For Mixtures of Restricted Two-Component Normals
#'
#' @param pi A vector of mixing probabilities.
#' @param p A vector of mixture weights as in [drtcnorm].
#' @param mu A vector of separation parameters.
#' @param sig1 A vector of standard deviations for the first component.
#' @param sig2 A vector of standard deviations for the second component.
#'
#' @return A variance scaling term.
#' @seealso [rtcnorm]
#' @keywords internal
var_e_fn <- function(pi, p, mu, sig1, sig2) {
  c0 <- p ^ 2 + (1 - p) ^ 2
  c1 <- (1 - p) / c0
  c2 <- -p / c0

  mu1 <- c1 * mu
  mu2 <- c2 * mu

  y <- p * (mu1 ^ 2 + sig1 ^ 2) + (1 - p) * (mu2 ^ 2 + sig2 ^ 2)
  y <- sum(pi * y)
  return(y)
}
