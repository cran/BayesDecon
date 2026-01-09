#' Generate a Proposal Variance Matrix for B-Spline Parameters
#'
#' @param thetas Spine parameters.
#' @param xs Data values.
#' @param mis Subject indicators.
#' @param us Residuals.
#' @param s2t Variance.
#' @param K.t Spline cardinality.
#' @param P.t Penalty matrix.
#' @param knots.t Knot points.
#' @param n Number of grid points.
#'
#' @return A proposal matrix.
#' @keywords internal
prop.sig.thetas.fn <- function(thetas, xs, mis, us, s2t, K.t, P.t, knots.t, n) {
  vars <- bspline_basis(xs, min(knots.t), max(knots.t), K.t) %*% exp(thetas)

  B.basis.components <- matrix(0, nrow = K.t + 1, ncol = n)
  for (kk in 1:(K.t + 1)) {
    thetas.new <- rep(0, K.t+1)
    thetas.new[kk] <- 1
    B.basis.components[kk, ] <- bspline_basis(xs, min(knots.t), max(knots.t),
                                              K.t) %*% thetas.new
  }

  prop.sig.thetas <- matrix(0, nrow = K.t + 1, ncol = K.t + 1)
  for (kk in 1:(K.t + 1)) {
    for (ll in kk:(K.t + 1)) {
      if (kk == ll) {
        for (ii in 1:n) {
          for (jj in (sum(mis[1:(ii - 1)]) + 1):sum(mis[1:ii])) {
            prop.sig.thetas[kk, ll] <- prop.sig.thetas[kk, ll] +
              ((us[jj] ^ 2) / vars[ii] - 1 / 2) *
              (B.basis.components[kk, ii] * B.basis.components[ll, ii]) *
              exp(thetas[kk] + thetas[ll]) / (vars[ii] ^ 2) +
              (1 - (us[jj] ^ 2) / vars[ii]) * B.basis.components[kk, ii] *
              exp(thetas[kk]) / (2 * vars[ii])
          }
        }
      }
      else {
        for (ii in 1:n) {
          for (jj in (sum(mis[1:(ii - 1)]) + 1):sum(mis[1:ii])) {
            prop.sig.thetas[kk, ll] <- prop.sig.thetas[kk, ll] +
              ((us[jj] ^ 2) / vars[ii] - 1 / 2) *
              (B.basis.components[kk, ii] * B.basis.components[ll, ii]) *
              exp(thetas[kk] + thetas[ll]) / (vars[ii] ^ 2)
          }
        }
      }
    }
  }
  for (kk in 2:(K.t + 1)) {
    for (ll in 1:(kk - 1)) {
      prop.sig.thetas[kk, ll] <- prop.sig.thetas[ll, kk]
    }
  }
  prop.sig.thetas <- prop.sig.thetas + P.t / s2t
  prop.sig.thetas <- round(solve(prop.sig.thetas), 4)

  return(prop.sig.thetas)
}




#' Generate a Proposal Variance Matrix for B-Spline Parameters for episodic density
#'
#' @param thetas Spine parameters.
#' @param xs Data values.
#' @param s2t Variance.
#' @param K.t Spline cardinality.
#' @param P.t Penalty matrix.
#' @param knots.t Knot points.
#' @param n Number of grid points.
#'
#' @return A proposal matrix.
#' @keywords internal
prop.sig.thetas.x.dens.fn <- function(thetas, xs, s2t, K.t, P.t, knots.t, n) {
  vars <- bspline_basis(xs, min(knots.t), max(knots.t), K.t) %*% exp(thetas)
  delta = knots.t[2] - knots.t[1]
  area_theta <- (delta / 6) * c(1,5,rep(6,K.t-3),5,1)
  vars_ <- sum(area_theta * exp(thetas))

  B.basis.components <- matrix(0, nrow = K.t + 1, ncol = n)
  for (kk in 1:(K.t + 1)) {
    thetas.new <- rep(0, K.t+1)
    thetas.new[kk] <- 1
    B.basis.components[kk, ] <- bspline_basis(xs, min(knots.t), max(knots.t),
                                              K.t) %*% thetas.new
  }

  prop.sig.thetas <- matrix(0, nrow = K.t + 1, ncol = K.t + 1)
  for (kk in 1:(K.t + 1)) {
    for (ll in kk:(K.t + 1)) {
      if (kk == ll) {
        for (ii in 1:n) {
          prop.sig.thetas[kk, ll] <- prop.sig.thetas[kk, ll] +
            (B.basis.components[kk, ii] * B.basis.components[ll, ii]) *
            exp(thetas[kk] + thetas[ll]) / (vars[ii] ^ 2) -
            B.basis.components[kk, ii] * exp(thetas[kk]) / vars[ii] -
            (area_theta[kk] * area_theta[ll]) *
            exp(thetas[kk] + thetas[ll]) / (vars_ ^ 2) +
            area_theta[kk] * exp(thetas[kk]) / vars_
        }
      }
      else {
        for (ii in 1:n) {
          prop.sig.thetas[kk, ll] <- prop.sig.thetas[kk, ll] +
            (B.basis.components[kk, ii] * B.basis.components[ll, ii]) *
            exp(thetas[kk] + thetas[ll]) / (vars[ii] ^ 2) -
            (area_theta[kk] * area_theta[ll]) *
            exp(thetas[kk] + thetas[ll]) / (vars_ ^ 2)
        }
      }
    }
  }
  for (kk in 2:(K.t + 1)) {
    for (ll in 1:(kk - 1)) {
      prop.sig.thetas[kk, ll] <- prop.sig.thetas[ll, kk]
    }
  }
  prop.sig.thetas <- prop.sig.thetas + P.t / s2t
  prop.sig.thetas <- prop.sig.thetas + diag(rep(1,K.t + 1))
  prop.sig.thetas <- round(solve(prop.sig.thetas), 4)
  return(prop.sig.thetas)
}
