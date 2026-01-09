#' Computing mean of MCMC samples
#' @param x A `bdeconv_result` object.
#' @param ... Unused.
#'
#' @return A `bdeconv_result` object with one sample constituting the mean
#  result of `x`.
#'
#' @details Where a direct point-wise mean is not appropriate (i.e. for mixture
#'  density results), mixtures are expanded to included every component from
#'  every sampled density with weights divided by `n`, the number of samples
#'  from which the mean is being taken.
#' @exportS3Method
mean.bdeconv_result <- function(x, ...) {
  obj <- x
  len <- length(obj)
  if (len == 1) {
    return(obj)
  } else {
    NextMethod()
  }
}

#' @exportS3Method
mean.zi_res <- function(x, ...) {
  obj <- x
  len <- length(obj)

  obj$w <- matrix(apply(obj$w, 2, mean), nrow = 1)

  obj$x <- matrix(apply(obj$x, 2, mean), nrow = 1)
  obj$theta_x <- matrix(t(apply(obj$theta_x, 2, mean)), nrow = 1)

  obj$x_e <- matrix(apply(obj$x_e, 2, mean), nrow = 1)
  obj$theta_e <- matrix(t(apply(obj$theta_e, 2, mean)), nrow = 1)

  obj$x_nz <- matrix(apply(obj$x_nz, 2, mean), nrow = 1)

  obj$theta_v <- matrix(t(apply(obj$theta_v, 2, mean)), nrow = 1)
  obj$var_e <- mean(obj$var_e)

  obj$u <- matrix(apply(obj$u, 2, mean), nrow = 1)
  K <- sum(obj$k_u)
  obj$k_u <- K
  obj$z_u <- NULL
  pi_u <- as.numeric(t(obj$pi_u))
  p_u <- as.numeric(t(obj$p_u))
  mu_u <- as.numeric(t(obj$mu_u))
  sig1_u <- as.numeric(t(obj$sig1_u))
  sig2_u <- as.numeric(t(obj$sig2_u))
  srt <- sort.int(pi_u, decreasing = TRUE, index.return = TRUE)
  idx <- srt$ix
  obj$pi_u <- srt$x / sum(srt$x)
  obj$p_u <- p_u[idx]
  obj$mu_u <- mu_u[idx]
  obj$sig1_u <- sig1_u[idx]
  obj$sig2_u <- sig2_u[idx]

  return(obj)
}


#' @exportS3Method
mean.nc_res <- function(x, ...) {
  obj <- x
  len <- length(obj)

  obj$x <- apply(obj$x, 2, mean)
  obj$z_x <- NULL
  pi_x <- as.numeric(t(obj$pi_x))
  mu_x <- as.numeric(t(obj$mu_x))
  sig_x <- as.numeric(t(obj$sig_x))
  srt <- sort.int(pi_x, decreasing = TRUE, index.return = TRUE)
  idx <- srt$ix
  obj$pi_x <- pi_x[idx] / len
  obj$mu_x <- mu_x[idx]
  obj$sig_x <- sig_x[idx]

  obj$theta_v <- t(apply(obj$theta_v, 2, mean))
  obj$var_e <- mean(obj$var_e)

  obj$u <- apply(obj$u, 2, mean)
  K <- sum(obj$k_u)
  KK <- c(1, cumsum(obj$k_u))
  obj$k_u <- K
  obj$z_u <- NULL
  pi_u <- numeric(K)
  p_u <- numeric(K)
  mu_u <- numeric(K)
  sig1_u <- numeric(K)
  sig2_u <- numeric(K)
  for (i in 1:length(obj)) {
    k_span <- KK[i]:KK[i + 1]
    pi_u[k_span] <- obj$pi_u[i, k_span]
    p_u[k_span] <- obj$p_u[i, k_span]
    mu_u[k_span] <- obj$mu_u[i, k_span]
    sig1_u[k_span] <- obj$sig1_u[i, k_span]
    sig2_u[k_span] <- obj$sig2_u[i, k_span]
  }
  srt <- sort.int(pi_u, decreasing = TRUE, index.return = TRUE)
  idx <- srt$ix
  obj$pi_u <- srt$x / sum(srt$x)
  obj$p_u <- p_u[idx]
  obj$mu_u <- mu_u[idx]
  obj$sig1_u <- sig1_u[idx]
  obj$sig2_u <- sig2_u[idx]

  return(obj)
}

#' @exportS3Method
mean.mv_res <- function(x, ...) {
  obj <- x
  len <- length(obj)
  d <- obj$p + obj$q

  obj$w <- apply(obj$w, 2:3, mean)
  dim(obj$w) <- c(1, dim(obj$w))

  obj$x <- apply(obj$x, 2:3, mean)
  dim(obj$x) <- c(1, dim(obj$x))
  obj$theta_x <- apply(obj$theta_x, 2:3, mean)
  dim(obj$theta_x) <- c(1, dim(obj$theta_x))

  obj$x_e <- apply(obj$x_e, 2:3, mean)
  dim(obj$x_e) <- c(1, dim(obj$x_e))
  obj$theta_e <- apply(obj$theta_e, 2:3, mean)
  dim(obj$theta_e) <- c(1, dim(obj$theta_e))

  obj$x_nz <- apply(obj$x_nz, 2:3, mean)
  dim(obj$x_nz) <- c(1, dim(obj$x_nz))

  obj$theta_v <- apply(obj$theta_v, 2:3, mean)
  dim(obj$theta_v) <- c(1, dim(obj$theta_v))
  obj$var_e <- matrix(apply(obj$var_e, 1, mean), nrow = 1)

  obj$z_x <- NULL
  pi_x  <- array(dim = c(1, obj$p, obj$control$K_x * len))
  mu_x  <- array(dim = c(1, obj$p, obj$control$K_x * len))
  sig_x <- array(dim = c(1, obj$p, obj$control$K_x * len))
  if (obj$p > 0) {
    for (pp in 1:obj$p) {
      pi_x[, pp, ]  <- as.numeric(t(obj$pi_x[, pp, ]))
      mu_x[, pp, ]  <- as.numeric(t(obj$mu_x[, pp, ]))
      sig_x[, pp, ] <- as.numeric(t(obj$sig_x[, pp, ]))
      srt <- sort.int(pi_x[, pp, ], decreasing = TRUE, index.return = TRUE)
      idx <- srt$ix
      pi_x[, pp, ]  <- pi_x[, pp, idx] / len
      mu_x[, pp, ]  <- mu_x[, pp, idx]
      sig_x[, pp, ] <- sig_x[, pp, idx]
    }
  }
  obj$pi_x  <- pi_x
  obj$mu_x  <- mu_x
  obj$sig_x <- sig_x

  obj$u <- apply(obj$u, 2:3, mean)
  dim(obj$u) <- c(1, dim(obj$u))
  K <- apply(obj$k_u, 2, sum)
  obj$k_u <- matrix(K, nrow = 1)
  obj$z_u <- NULL
  pi_u <- t(sapply(1:d, function(dd) as.numeric(t(obj$pi_u[, dd, ]))))
  p_u <- as.numeric(t(obj$p_u))
  mu_u <- as.numeric(t(obj$mu_u))
  sig1_u <- as.numeric(t(obj$sig1_u))
  sig2_u <- as.numeric(t(obj$sig2_u))
  srt <- sort.int(apply(pi_u, 2, sum), decreasing = TRUE, index.return = TRUE)
  idx <- srt$ix
  obj$pi_u <- pi_u[, idx] / apply(pi_u, 1, sum)
  dim(obj$pi_u) <- c(1, dim(obj$pi_u))
  obj$p_u <- matrix(p_u[idx], nrow = 1)
  obj$mu_u <- matrix(mu_u[idx], nrow = 1)
  obj$sig1_u <- matrix(sig1_u[idx], nrow = 1)
  obj$sig2_u <- matrix(sig2_u[idx], nrow = 1)

  obj$corr_x <- apply(obj$corr_x, 2:3, mean)
  dim(obj$corr_x) <- c(1, dim(obj$corr_x))
  obj$corr_e <- apply(obj$corr_e, 2:3, mean)
  dim(obj$corr_e) <- c(1, dim(obj$corr_e))

  return(obj)
}
