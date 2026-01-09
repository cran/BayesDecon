#' The Truncated Normal Distribution
#'
#' Density, distribution, quantile, and random variate generation functions for
#'  the Truncated Normal distribution, with `mu` and `sig` equal to the mean and
#'  standard deviation before truncation and truncated to `[lwr, upr]`. Wrappers
#'  for functions from `msm`.
#'
#' @aliases tnorm
#' @aliases dtnorm
#' @aliases ptnorm
#' @aliases qtnorm
#' @aliases rtnorm
#'
#' @name TruncNormal
#'
#' @param x,q A vector of quantiles.
#' @param p A vector of probabilities.
#' @param n A number of observations.
#' @param mu A real-valued mean.
#' @param sig A strictly positive standard deviation.
#' @param lwr A lower bound.
#' @param upr An upper bound.
#' @param log,log.p A logical value indication whether a probability `p` should be
#'  given as `log(p)`
#' @param lower.tail A logical value indicating whether to take lower or upper
#'  tail probabilities.
#'
#' @return `dtnorm` gives the density, `ptnorm` gives the distribution function,
#'  `qtnorm` gives the quantile function, and `rtnorm` generates random
#'  variates.
#'  The length of the result is determined by the length of `x` for `dtnorm`.

#' @rdname tnorm
#' @keywords internal
dtnorm <- function(x, mu, sig, lwr, upr, log = FALSE) {
  msm::dtnorm(x, mu, sig, lwr, upr, log)
}

#' @rdname tnorm
#' @keywords internal
ptnorm <- function(q, mu, sig, lwr, upr, lower.tail = TRUE, log.p = FALSE) {
  msm::ptnorm(q, mu, sig, lwr, upr, lower.tail, log.p)
}

#' @rdname tnorm
#' @keywords internal
qtnorm <- function(p, mu, sig, lwr, upr, lower.tail = TRUE, log.p = FALSE) {
  msm::qtnorm(p, mu, sig, lwr, upr, lower.tail, log.p)
}

#' @rdname tnorm
#' @keywords internal
rtnorm <- function (n, mu, sig, lwr = -Inf, upr = Inf) {
  msm::rtnorm(n, mu, sig, lwr, upr)
}
