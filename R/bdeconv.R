#' Obtaining MCMC samples from deconvoluted density
#'
#' @param obj A data frame where the first column is subject labels.
#' @param lwr Unused.
#' @param upr The ratio of an upper bound for the deconvolved estimates to the upper bound of the surrogates. It can be either a
#'  scalar or a vector with a length that matches the number of numeric columns
#'   in `obj`.
#' @param wgts Unused.
#' @param regular_var Specifies which variables to consider as regular; all others are treated as zero inflated.
#' @param res_type The type of result returned from the MCMC sampler.
#' @param err_type The shape of the error distribution used for deconvolution.
#' @param var_type whether to use homogeneous or heterogeneous error distribution.
#' @param na_rm A logical value indicating whether to ignore rows which contain
#'  `NA` values in numeric columns.
#' @param control A list of arguments specifying prior hyper-parameters and parameters of the MCMC sampler.
#' @param progress Whether a progress bar for MCMC will be shown or not.  Default is FALSE.
#' @param parallel For multivariate data. Whether the pre-processing of each variable will be parallelized or not. Default is FALSE.
#' @param core_num If parallelized, number of cores to be used.
#' @return A `bdeconv_result` object, containing all the posterior samples, suitable for further analysis and visualization.
#'
#' @details All hyper-parameters of the model and variables necessary for the MCMC sampler can be specified through a list via the control argument; some of the most important ones are listed.
#'
#' - `nsamp`: Number of finally retained MCMC samples after burn-in and thinning. Default is 1000.
#' - `nburn`: Burn-in period. Default is 2000.
#' - `nthin`: Thinning interval. Default is 5.
#' - `niter_uni`: For multivariate data, number of iteration for pre-processing of each variable. Default is 500.
#' - `K_t`: Number of knot points for bspline. Default is 10.
#' - `K_x`: Number of mixture components in the model for the density of X when the surrogates are strictly continuous. Not applicable when the surrogates are zero-inflated in which case the density is modeled using normalized bsplines. Default is 10.
#' - `K_epsilon`: Number of mixture components in the model for the density of scaled errors when applicable. Default is 10.
#' - `sigmasq_epsilon`: variance parameter for the normal prior of the scaled error mixture components' location parameter. Default is 4.
#' - `a_epsilon`: shape parameter for the inverse gamma prior of the variances of scaled error mixture components. Default is 3.
#' - `b_epsilon`: rate parameter for the inverse gamma prior of the variances of scaled error mixture components. Default is 2.
#' - `a_vartheta`: shape parameter for the inverse gamma prior of the variances of bspline coefficients corresponding to the variance function. Default is 100.
#' - `b_vartheta`: rate parameter for the inverse gamma prior of the variances of bspline coefficients corresponding to the variance function. Default is 1.
#' - `a_xi`: shape parameter for the inverse gamma prior of the variances of bspline coefficients corresponding to the densities of episodic components. Default is 100.
#' - `b_xi`: rate parameter for the inverse gamma prior of the variances of bspline coefficients corresponding to the densities of episodic components. Default is 10.
#' - `a_beta`: shape parameter for the inverse gamma prior of the variances of bspline coefficients corresponding to the probability of consumption function. Default is 100.
#' - `b_beta`: rate parameter for the inverse gamma prior of the variances of bspline coefficients corresponding to the probability of consumption function. Default is 1.
#'
#' @rdname bdeconv
#' @examples
#' data(BayesDecon_data_simulated)
#' set.seed(123)
#' ### Selecting first 300 rows for demonstration.
#' ### For best results, using the whole data is recommended.
#' result_uni <- bdeconv(BayesDecon_data_simulated[1:300,c(1,2)],
#'  upr = 1, res_type = "all",
#'  control = list(nsamp = 500, nburn = 0, nthin = 1),
#'  progress = TRUE)
#' plot(result_uni)
#'
#' \donttest{
#' set.seed(123)
#' result_mult <- bdeconv(BayesDecon_data_simulated, upr = 1, res_type = "all", progress = TRUE)
#' plot(result_mult)
#' plot(result_mult[,c(1,2)])
#' value = rbind(c(1,1,1),c(2,2,2),c(3,3,3),c(4,4,4))
#' colMeans(dist_x(result_mult, vals = value))
#' colMeans(dens_x(result_mult, vals = value))
#' print(result_mult)
#' }
#'
# @seealso [mean.bdeconv_result()]
#' @export
bdeconv <- function(obj, lwr = 0, upr = NULL, wgts = NULL,
                    regular_var = NULL,
                    res_type = c('final', 'mean', 'all'),
                    err_type = c('flex', 'norm', 'lapl'),
                    var_type = c('heteroscedastic', 'homoscedastic'),
                    na_rm = FALSE, control = list(), progress = FALSE, core_num = 1, parallel = FALSE) {
  UseMethod('bdeconv')
}


#' @export
bdeconv.data.frame <- function(obj, lwr = 0, upr = NULL, wgts = NULL,
                               regular_var = NULL,
                               res_type = c('final', 'mean', 'all'),
                               err_type = c('flex', 'norm', 'lapl'),
                               var_type = c('heteroscedastic', 'homoscedastic'),
                               na_rm = FALSE, control = list(), progress = FALSE, core_num = 1, parallel = FALSE) {
  if (!is.null(wgts)) {
    excl <- c(1, which(names(obj) == wgts))
  } else {
    excl <- 1
  }
  y <- obj[, -excl, drop = FALSE]
  y <- y[sapply(y, is.numeric)]
  y_org <- y
  if (ncol(y) == 0) {
    stop('`obj` must have at least one numeric column.')
  }
  if (is.null(upr)) {
    upr <- 0.5
  }
  if (length(upr) == 1) {
    upr <- rep(upr, ncol(y))
  }
  if (!is.null(lwr) & length(lwr) == 1) {
    lwr <- rep(lwr, ncol(y))
  }
  scl <- rep(1, ncol(y))
  loc <- rep(0, ncol(y))
  for (col in 1:ncol(y)) {
    ymax <- max(y[, col], na.rm = TRUE)
    ymin <- max(min(y[, col], na.rm = TRUE),0)
    if (ymin > 0){
      regular_var <- c(regular_var, colnames(y)[col])
    }
    if (TRUE) {
      scl[col] <- 10 / (ymax-ymin)
      loc[col] <- ymin
      y[, col] <- 10 * (y[, col]-ymin) / (ymax-ymin)
      upr[col] <- 10 * upr[col]
    } else {
      upr[col] <- ymax * upr[col]
    }

  }
  regular_var <- unique(regular_var)
  obj[, -excl] <- y
  if (!is.null(regular_var)){
    if (!any(regular_var %in% names(y))){
      regular_var <- NULL
      warning("argument regular_var doesn't matches with column names of object, assigning regular_var to NULL")
    } else{
      for (varname in regular_var){
        if(varname %in% names(y)){
        if (any(y[,varname]==0)){
          replace_value=min(y[y[,varname]>0,varname])/2
          y[,varname]<-replace(y[,varname],y[,varname]==0,replace_value)
        }
        }
      }
    }
  }
  obj[, -excl] <- y
  if (ncol(y) == 1) {
    if (any(y == 0)) {
      class(obj) <- c('unizi', class(obj))
    } else {
      class(obj) <- c('uninc', class(obj))
    }
  } else {
    idx_zi <- which(apply(y, 2, function(c) any(c == 0)))
    idx_nc <- which(apply(y, 2, function(c) !any(c == 0)))
    q <- length(idx_zi)
    p <- length(idx_nc)
    d <- p + q
    y <- y[, c(idx_zi, idx_nc)]
    y_org <- y_org[, c(idx_zi, idx_nc)]
    upr <- upr[c(idx_zi, idx_nc)]
    lwr <- lwr[c(idx_zi, idx_nc)]
    class(obj) <- c('mvtzi', class(obj))
  }

  ## checking whether parallel computing is possible or not
  if (parallel) {
    # Check if required packages are installed and loadable
    pkgs_needed <- c("parallel", "foreach", "doParallel")
    pkgs_available <- vapply(pkgs_needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))

    if (!all(pkgs_available)) {
      warning("Some packages for parallel processing are missing. Setting parallel = FALSE.")
      parallel <- FALSE
    }
  }

  res <- bdeconv(obj, lwr, upr, wgts,regular_var, res_type, err_type, var_type, na_rm,
                 control, progress, core_num, parallel)
  res <- bd_descale(res, 1 / scl, loc)
  res$scl <- scl
  res$names <- names(y)
  res$call <- match.call()
  res$y <- as.matrix(y_org)
  return(res)
}


#' @export
bdeconv.unizi <- function(obj, lwr = 0, upr = NULL, wgts = NULL,
                          regular_var = NULL,
                          res_type = c('final', 'mean', 'all'),
                          err_type = c('flex', 'norm', 'lapl'),
                          var_type = c('heteroscedastic', 'homoscedastic'),
                          na_rm = FALSE, control = list(), progress = FALSE, core_num = 1, parallel = FALSE) {
  # Unpack object and validate data
  if (!is.null(wgts)) {
    excl <- c(1, which(names(obj) == wgts))
  } else {
    excl <- 1
  }
  m <- table(obj[, 1])
  w <- obj[, -excl]
  n <- length(m)
  N <- sum(m)
  w.na <- is.na(w)
  if (any(w.na)) {
    if (na_rm) {
      N <- length(w <- w[!w.na])
    }
    else {
      stop('`w` contains missing values.')
    }
  }
  N <- as.integer(N)
  if (is.na(N)) {
    stop(gettextf('Invalid value of %s', 'length(w).'), domain = NA)
  }
  if (any(!is.finite(w))) {
    stop('`w` contains non-finite values.')
  }
  # Get model customizations
  err_type <- match.arg(err_type)
  err_type_c <- which(c('flex', 'norm', 'lapl') == err_type) - 1
  var_type <- match.arg(var_type)
  var_type_c <- var_type == 'heteroscedastic'
  # Bind unbound variables
  nsamp <- nburn <- nthin <- niter_u <- NULL
  x_grd_res <- var_grd_res <- e_grd_res <- NULL
  K_t <- K_epsilon <- a_vartheta <- b_vartheta <- sigmasq_vartheta <- NULL
  a_beta <- b_beta <- sigmasq_beta <- NULL
  a_xi <- b_xi <- sigmasq_xi <- NULL
  alpha_epsilon <- a_epsilon <- b_epsilon <- sigmasq_epsilon <- NULL
  sigma00_theta_e <- sig_tune_theta_1 <- sig_tune_theta_2 <- NULL
  # Default control parameters (hyperparameters, sampler parameters, etc.)
  control_defs <- list(nsamp = 1000, nburn = 2000, nthin = 5, niter_u = 10,
                       x_grd_res = 500, var_grd_res = 100, e_grd_res = 100,
                       K_t = 10, K_epsilon = 10,
                       a_vartheta = 100, b_vartheta = 1, sigmasq_vartheta = 0.01,
                       a_beta = 100, b_beta = 1, sigmasq_beta = 0.01,
                       a_xi = 100, b_xi = 10,
                       sigmasq_xi = 0.1, alpha_epsilon = .1, a_epsilon = 3, b_epsilon = 2,
                       sigmasq_epsilon = 4, ## sigmasq_epsilon is actually the variance not sd
                       sigma00_theta_e = 0.1,
                       sig_tune_theta_1 = 0, sig_tune_theta_2 = 0.1)
  for (nm in names(control_defs)) {
    if (is.null(control[[nm]])) {
      control[[nm]] <- control_defs[[nm]]
    }
    assign(nm, control[[nm]])
  }
  # Validate control parameters
  validate_sclr_params(nsamp, nburn, nthin, niter_u,
                       x_grd_res, var_grd_res, e_grd_res)
  # Initialize `x` and `u`
  idx <- match(obj[, 1], sort(unique(obj[, 1])))
  w_bar <- tapply(w, idx, mean)
  w_var <- tapply(w, idx, stats::var)
  x <- as.numeric(w_bar)
  x[x <= lwr] <- lwr + 0.1
  x[x >= upr] <- upr - 0.1
  u <- w - rep(x, times = m)
  # Initialize thetas for variance function
  P2_t <- P2.mat(K_t + 1)
  knots_t <- seq(lwr, upr, length = K_t)
  optim_results <- stats::optim(rep(1, K_t + 1), fr, NULL,
                                x, m, knots_t, P2_t, sigmasq_vartheta, u,
                                method = 'BFGS')
  theta <- optim_results$par
  prop_sig_theta <- corpcor::make.positive.definite(
    prop.sig.thetas.fn(theta, x, m, u, sigmasq_vartheta, K_t, P2_t, knots_t, n)
  )
  if (!var_type_c) {
    theta <- rep(0, K_t + 1)
  }
  x_grd <- seq(lwr, upr, length = x_grd_res)
  kde_dens_est <- ks::kde(x, eval.points = x_grd)$estimate
  optim_results <- stats::optim(rep(1, K_t + 1), fr.dens, NULL,
                                x_grd, knots_t, P2_t, sigmasq_xi, kde_dens_est,
                                K_t,
                                method = 'BFGS')
  theta_x_dens <- optim_results$par
  prop_sig_theta_x_dens <- corpcor::make.positive.definite(
    prop.sig.thetas.x.dens.fn(theta_x_dens, x, sigmasq_xi, K_t, P2_t, knots_t, n)
  )
  idx_nonzero <- which(w != 0)
  ybin <- numeric(N)
  ybin[idx_nonzero] <- 1
  probs_emp <- tapply(ybin, idx, sum) / m
  optim_results <- stats::optim(numeric(K_t + 1), fr.prob.nonzero, NULL,
                                w_bar, knots_t, P2_t, sigmasq_beta, probs_emp,
                                method = 'BFGS')
  mu0_theta_e <- cumsum(exp(optim_results$par))
  sigma0_theta_e <- sigma00_theta_e * diag(K_t + 1)
  # Conduct MCMC
  if(progress){
    cat(sprintf("\r Univariate MCMC chain started \n"))
    utils::flush.console()
  }
  res <- mcmcuniepi_cpp(nsamp, nburn, nthin, niter_u,
                        w, x, u, m, idx - 1,
                        lwr, upr,
                        var_type_c, K_t, knots_t, theta, P2_t,
                        x_grd_res, e_grd_res, var_grd_res,
                        theta_x_dens,
                        err_type_c, K_epsilon,
                        a_vartheta, b_vartheta, sigmasq_vartheta,
                        a_beta, b_beta, sigmasq_beta,
                        a_xi, b_xi, sigmasq_xi,
                        alpha_epsilon, a_epsilon, b_epsilon, sigmasq_epsilon,  ## sigmasq_epsilon is actually the variance not sd
                        mu0_theta_e, sigma0_theta_e,
                        prop_sig_theta_x_dens, prop_sig_theta, sig_tune_theta_1, sig_tune_theta_2,progress)
  # Process results
  res_type <- match.arg(res_type)
  res$q <- 1
  res$p <- 0
  res$w_var <- w_var
  res$lwr <- lwr
  res$upr <- upr
  res$control <- control
  res$res_type <- res_type
  res$err_type <- err_type
  return(bdeconvzi_result(res, res_type))
}

#' @export
bdeconv.uninc <- function(obj, lwr = 0, upr = NULL, wgts = NULL,
                          regular_var = NULL,
                          res_type = c('final', 'mean', 'all'),
                          err_type = c('flex', 'norm', 'lapl'),
                          var_type = c('heteroscedastic', 'homoscedastic'),
                          na_rm = FALSE, control = list(), progress = FALSE, core_num = 1, parallel = FALSE) {
  if (!is.null(wgts)) {
    excl <- c(1, which(names(obj) == wgts))
    wgts <- obj[, wgts]
  } else {
    excl <- 1
    wgts <- rep(1, nrow(obj))
  }
  # Unpack object and validate data
  m <- table(obj[, 1])
  w <- obj[, -excl]
  n <- length(m)
  N <- sum(m)
  w.na <- is.na(w)
  if (any(w.na)) {
    if (na_rm) {
      N <- length(w <- w[!w.na])
    }
    else stop('`w` contains missing values.')
  }
  N <- as.integer(N)
  if (is.na(N)) {
    stop(gettextf('Invalid value of %s', 'length(w).'), domain = NA)
  }
  if (any(!is.finite(w))) {
    stop('`w` contains non-finite values.')
  }
  # Get model customizations
  err_type <- match.arg(err_type)
  err_type_c <- which(c('flex', 'norm', 'lapl') == err_type) - 1
  var_type <- match.arg(var_type)
  var_type_c <- var_type == 'heteroscedastic'
  # Bind unbound variables
  nsamp <- nburn <- nthin <- niter_u <- NULL
  x_grd_res <- var_grd_res <- e_grd_res <- NULL
  K_t <- K_epsilon <- alpha_x <- a_sigmasq_x <- b_sigmasq_x <- K_x <- NULL
  a_vartheta <- b_vartheta <- sigmasq_vartheta <- NULL
  alpha_epsilon <- a_epsilon <- b_epsilon <- sigmasq_epsilon <- NULL
  sigma00_theta_e <- sig_tune_theta_1 <- sig_tune_theta_2 <- NULL
  # Default control parameters (hyperparameters, sampler parameters, etc.)
  control_defs <- list(nsamp = 1000, nburn = 2000, nthin = 5, niter_u = 10,
                       x_grd_res = 500, var_grd_res = 200, e_grd_res = 100,
                       K_t = 10, K_x = 10, K_epsilon = 10,
                       alpha_x = .1, a_sigmasq_x = 1, b_sigmasq_x = 1,
                       a_vartheta = 100, b_vartheta = 1, sigmasq_vartheta = 0.01,
                       alpha_epsilon = .1, a_epsilon = 3, b_epsilon = 2,
                       sigmasq_epsilon = 4,
                       sigma00_theta_e = 0.1,
                       sig_tune_theta_1 = 0, sig_tune_theta_2 = 0.1)
  for (nm in names(control_defs)) {
    if (is.null(control[[nm]])) {
      control[[nm]] <- control_defs[[nm]]
    }
    assign(nm, control[[nm]])
  }
  # Validate control parameters
  validate_sclr_params(nsamp, nburn, nthin, niter_u,
                       x_grd_res, var_grd_res, e_grd_res)
  # Initialize `x` and `u`
  idx <- match(obj[, 1], sort(unique(obj[, 1])))
  w_bar <- tapply(w, idx, mean)
  w_var <- tapply(w, idx, stats::var)
  x <- as.numeric(w_bar)
  x[x <= lwr] <- lwr + 0.1
  x[x >= upr] <- upr - 0.1
  u <- w - rep(x, times = m)
  # Initialize `x` clustering parmameters
  pi_x <- rep(1 / K_x, K_x)
  clusters_xs <- suppressWarnings(stats::kmeans(x, K_x, iter.max =50))
  mu_xs <- clusters_xs$centers
  sigmasq_xs <- rep(stats::var(x) / 5, K_x)
  z_xs <- clusters_xs$cluster
  ## Initialize `x` normal prior paramaters
  mu0_x <- mean(x)
  sigmasq0_x <- stats::var(x)
  # Initialize thetas for variance function
  P_t <- P2.mat(K_t + 1)
  knots_t <- seq(lwr, upr, length = K_t)
  optim_results <- stats::optim(rep(1,K_t + 1), fr, NULL,
                                x, m, knots_t, P_t, sigmasq_vartheta, u,
                                method = 'BFGS')
  theta <- optim_results$par
  prop_sig_theta <- corpcor::make.positive.definite(
    prop.sig.thetas.fn(theta, x, m, u, sigmasq_vartheta, K_t, P_t, knots_t, n)
  )
  if (!var_type_c) {
    theta <- rep(0, K_t + 1)
  }
  # Conduct MCMC
  if(progress){
    cat(sprintf("\r Univariate MCMC chain started \n"))
    utils::flush.console()
  }
  res <- mcmcunireg_cpp(nsamp, nburn, nthin, niter_u,
                        w, x, u, m, idx - 1, lwr, upr, wgts,
                        var_type_c, K_t,
                        knots_t, theta, P_t,
                        x_grd_res, e_grd_res, var_grd_res,
                        K_x, pi_x, z_xs, mu_xs, sigmasq_xs,
                        err_type_c, K_epsilon, alpha_x, alpha_epsilon,
                        mu0_x, sigmasq0_x, a_sigmasq_x, b_sigmasq_x,
                        a_epsilon, b_epsilon, sigmasq_epsilon,
                        sigmasq_vartheta, a_vartheta, b_vartheta,
                        prop_sig_theta, sig_tune_theta_1, sig_tune_theta_2,
                        FALSE,progress)
  # Process results
  res_type <- match.arg(res_type)
  res$w <- w
  res$q <- 0
  res$p <- 1
  res$w_var <- w_var
  res$lwr <- lwr
  res$upr <- upr
  res$control <- control
  res$res_type <- res_type
  res$err_type <- err_type
  return(bdeconvnc_result(res, res_type))
}

#' @export
bdeconv.mvtzi <- function(obj, lwr = 0, upr = NULL, wgts = NULL,
                          regular_var = NULL,
                          res_type = c('final', 'mean', 'all'),
                          err_type = c('flex', 'norm', 'lapl'),
                          var_type = c('heteroscedastic', 'homoscedastic'),
                          na_rm = FALSE, control = list(), progress = FALSE, core_num =1, parallel = FALSE) {
  # Unpack object and validate data
  inds <- obj[, 1]
  m <- table(obj[, 1])
  w <- obj[, -1]
  n <- length(m)
  N <- sum(m)
  w.na <- is.na(w)
  if (any(w.na)) {
    if (na_rm) {
      N <- length(w <- stats::na.omit(w))
    }
    else stop('`w` contains missing values.')
  }
  N <- as.integer(N)
  if (is.na(N)) {
    stop(gettextf('Invalid value of %s', 'length(w).'), domain = NA)
  }
  if (any(!apply(w, 2, is.finite))) {
    stop('`w` contains non-finite values.')
  }
  # Find zero-inflated columns
  idx_zi <- which(apply(w, 2, function(c) any(c == 0)))
  idx_nc <- which(apply(w, 2, function(c) !any(c == 0)))
  q <- length(idx_zi)
  p <- length(idx_nc)
  d <- p + q
  w <- w[, c(idx_zi, idx_nc)]
  upr <- upr[c(idx_zi, idx_nc)]
  lwr <- lwr[c(idx_zi, idx_nc)]
  # Get model customizations
  err_type <- match.arg(err_type)
  err_type_c <- which(c('flex', 'norm', 'lapl') == err_type) - 1
  var_type <- match.arg(var_type)
  var_type_c <- var_type == 'heteroscedastic'
  # Bind unbound variables
  nsamp <- nburn <- nthin <- niter_u <- niter_uni <- NULL
  x_grd_res <- var_grd_res <- e_grd_res <- NULL
  K_t <- K_epsilon <- alpha_x <- a_sigmasq_x <- b_sigmasq_x <- K_x <- NULL
  a_vartheta <- b_vartheta <- sigmasq_vartheta <- a_beta <- b_beta <- sigmasq_beta <- NULL
  a_xi <- b_xi <- sigmasq_xi <- NULL
  alpha_epsilon <- a_epsilon <- b_epsilon <- sigmasq_epsilon <- NULL
  sigma00_theta_e <- mu_theta_e <- sig_tune_theta_1 <- sig_tune_theta_2 <- NULL
  # Default control parameters (hyperparameters, sampler parameters, etc.)
  control_defs <- list(nsamp = 1000, nburn = 2000, nthin = 5, niter_u = 10,
                       niter_uni = 500,
                       x_grd_res = 500, var_grd_res = 200, e_grd_res = 100,
                       K_t = 10, K_epsilon = 10,
                       alpha_x = 0.1, a_sigmasq_x = 1, b_sigmasq_x = 1,
                       K_x = 10,
                       a_vartheta = 100, b_vartheta = 1, sigmasq_vartheta = 0.01,
                       a_beta =  rep(100, q), b_beta = rep(1, q),
                       sigmasq_beta = rep(0.01, q),
                       a_xi = rep(100, q), b_xi =  rep(10, q),
                       sigmasq_xi = rep(0.1, q),
                       alpha_epsilon = .1, a_epsilon = 3, b_epsilon = 2,
                       sigmasq_epsilon = 4,
                       sigma00_theta_e = 0.1,
                       mu0_theta_e = matrix(0, q, 10 + 1),
                       sig_tune_theta_1 = 0, sig_tune_theta_2 = .1)
  for (nm in names(control_defs)) {
    if (is.null(control[[nm]])) {
      control[[nm]] <- control_defs[[nm]]
    } else{
      if (nm %in% c('a_beta', 'b_beta', 'sigmasq_beta', 'a_xi', 'b_xi', 'sigmasq_xi') & length(control[[nm]]) == 1){
        control[[nm]] <- rep(control[[nm]] , q)
      }
    }
    assign(nm, control[[nm]])
  }
  if (ncol(mu0_theta_e) != K_t + 1) {
    if (all(mu0_theta_e == 0)) {
      mu0_theta_e <- matrix(0, q, K_t + 1)
    } else {
      stop('`mu0_theta_e` must have `K_t` + 1 columns. Check your `control` ',
           'params.')
    }
  }

  w <- t(w)
  w_var <- apply(w, 1, function(ww) tapply(ww, inds, stats::var))
  x <- matrix(0, nrow = d, ncol = n)
  u <- matrix(0, nrow = d, ncol = N)
  knots <- matrix(0, nrow = d, ncol = K_t)
  delta <- numeric(d)
  theta_v <- matrix(0, nrow = d, ncol = K_t + 1)
  theta_x <- matrix(0, nrow = q, ncol = K_t + 1)
  z_u <- matrix(nrow = d, ncol = N)
  pi_u <- matrix(0, nrow = d, ncol = K_epsilon)
  params_u <- array(0, dim = c(d, K_epsilon, 4))
  ic <- w != 0
  inc <- 1 - ic
  x_e <- matrix(0, nrow = q, ncol = n)
  x_nz <- matrix(0, nrow = d, ncol = n)
  theta_e <- matrix(0, nrow = d, ncol = K_t + 1)
  mu0_theta_e <- matrix(0, nrow = q, K_t + 1)
  z_x <- matrix(0, nrow = p, ncol = n)
  pi_x <- matrix(0, nrow = p, ncol = K_x)
  params_x <- array(0, dim = c(p, K_x, 2))
  prop_sig_theta <- array(0, dim = c(d, K_t + 1, K_t + 1))
  prop_sig_theta_x_dens <- array(0, dim = c(q, K_t + 1, K_t + 1))
  P_t <- P2.mat(K_t + 1)
  # Episodic variable preprocessing
  if (parallel & q > 0) {
    # Detect cores and create cluster
    num_cores <- max(min(parallel::detectCores() - 1,core_num),1)
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    foreach_obj <- foreach::foreach(qq = 1:q, .packages = 'BayesDecon')
    results <-  foreach::`%dopar%`(foreach_obj, {
      obj_tmp <- obj[, c(1, (1 + qq))]
      class(obj_tmp)[[1]] <- 'unizi'

      control_uni <- control
      control_uni$nsamp <- niter_uni
      control_uni$nthin <- 1
      control_uni$nburn <- 0

      for (nm in c('a_beta', 'b_beta', 'sigmasq_beta', 'a_xi',
                   'b_xi', 'sigmasq_xi')) {
        control_uni[[nm]] <- control_uni[[nm]][qq]
      }

      control_uni$mu0_theta_e <- control_uni$mu0_theta_e[qq, ]

      # Progress off inside parallel loop for stability
      res <- bdeconv(obj_tmp, lwr = lwr[qq], upr = upr[qq],
                     res_type = 'final', err_type = err_type, var_type = var_type,
                     control = control_uni, progress = FALSE)

      list(
        w = res$w,
        x = res$x,
        u = res$u,
        knots = res$knots,
        delta = res$knots[2] - res$knots[1],
        theta_v = res$theta_v,
        theta_x = res$theta_x,
        z_u = res$z_u,
        pi_u = res$pi_u,
        params_u = cbind(res$p_u, res$mu_u, res$sig1_u, res$sig2_u),
        x_e = res$x_e,
        x_nz = res$x_nz,
        theta_e = res$theta_e,
        mu0_theta_e = res$theta_e
      )
    })
    parallel::stopCluster(cl)
    # Now update pre-existing matrices by rows
    for (qq in seq_along(results)) {
      w[qq, ] <- results[[qq]]$w
      x[qq, ] <- results[[qq]]$x
      u[qq, ] <- results[[qq]]$u
      knots[qq, ] <- results[[qq]]$knots
      delta[qq] <- results[[qq]]$delta
      theta_v[qq, ] <- results[[qq]]$theta_v
      theta_x[qq, ] <- results[[qq]]$theta_x
      z_u[qq, ] <- results[[qq]]$z_u
      pi_u[qq, 1:K_epsilon] <- results[[qq]]$pi_u
      params_u[qq, 1:K_epsilon, ] <- results[[qq]]$params_u
      x_e[qq, ] <- results[[qq]]$x_e
      x_nz[qq, ] <- results[[qq]]$x_nz
      theta_e[qq, ] <- results[[qq]]$theta_e
      mu0_theta_e[qq, ] <- results[[qq]]$mu0_theta_e
    }
  } else if (!parallel){
    for (qq in 1:q) {
      if (q == 0) {
        next
      }
      obj_tmp <- obj[, c(1, (1 + qq))]
      class(obj_tmp)[[1]] <- 'unizi'
      control_uni <- control
      control_uni$nsamp <- niter_uni
      control_uni$nthin <- 1
      control_uni$nburn <- 0
      for (nm in c('a_beta', 'b_beta', 'sigmasq_beta', 'a_xi',
                   'b_xi', 'sigmasq_xi')) {
        control_uni[[nm]] <- control_uni[[nm]][qq]
      }
      control_uni$mu0_theta_e <- control_uni$mu0_theta_e[qq, ]
      if (progress){
        msg <- sprintf("preprocessing variable %s", names(obj_tmp)[2])
        cat(sprintf("\r%-*s", 50, msg))  # left-align + pad with spaces
        utils::flush.console()
      }
      res <- bdeconv(obj_tmp, lwr = lwr[qq], upr = upr[qq],
                     res_type = 'final', err_type = err_type, var_type = var_type,
                     control = control_uni, progress = FALSE)

      w[qq, ] <- res$w
      x[qq, ] <- res$x
      u[qq, ] <- res$u
      knots[qq, ] <- res$knots
      delta[qq] <- res$knots[2] - res$knots[1]
      theta_v[qq, ] <- res$theta_v
      theta_x[qq, ] <- res$theta_x

      z_u[qq, ] <- res$z_u
      pi_u[qq, 1:K_epsilon] <- res$pi_u
      params_u[qq, 1:K_epsilon, ] <- cbind(res$p_u, res$mu_u, res$sig1_u,
                                         res$sig2_u)
      x_e[qq, ] <- res$x_e
      x_nz[qq, ] <- res$x_nz
      theta_e[qq, ] <- res$theta_e
      mu0_theta_e[qq, ] <- res$theta_e
    }
  }
  # Regular variable preprocessing
  if(parallel & p>0){
    # Detect cores and create cluster
    num_cores <- max(min(parallel::detectCores() - 1,core_num),1)
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    foreach_obj <- foreach::foreach(pp = (q + 1):(q + p), .packages = 'BayesDecon')
    results <- foreach::`%dopar%`(foreach_obj, {
      obj_tmp <- obj[, c(1, (1 + pp))]
      class(obj_tmp)[[1]] <- 'uninc'

      control_uni <- control
      control_uni$nsamp <- niter_uni
      control_uni$nthin <- 1
      control_uni$nburn <- 0

      res <- bdeconv(obj_tmp, lwr = lwr[pp], upr = upr[pp],
                     res_type = 'final', err_type = err_type, var_type = var_type,
                     control = control_uni, progress = FALSE)

      list(
        x = res$x,
        u = res$u,
        knots = res$knots,
        delta = res$knots[2] - res$knots[1],
        theta_v = res$theta_v,

        z_x = res$z_x,
        pi_x = res$pi_x,
        mu_x = res$mu_x,
        sig_x = res$sig_x,

        z_u = res$z_u,
        pi_u = res$pi_u,
        params_u = cbind(res$p_u, res$mu_u, res$sig1_u, res$sig2_u),

        x_nz = res$x
      )
    })

    # Stop cluster after parallel processing
    parallel::stopCluster(cl)

    # Update matrices row-wise after parallel loop
    for (idx in seq_along(results)) {
      pp <- (q + idx)
      x[pp, ] <- results[[idx]]$x
      u[pp, ] <- results[[idx]]$u
      knots[pp, ] <- results[[idx]]$knots
      delta[pp] <- results[[idx]]$delta
      theta_v[pp, ] <- results[[idx]]$theta_v

      z_x[pp - q, ] <- results[[idx]]$z_x
      pi_x[pp - q, 1:K_x] <- results[[idx]]$pi_x
      params_x[pp - q, 1:K_x, 1] <- results[[idx]]$mu_x
      params_x[pp - q, 1:K_x, 2] <- results[[idx]]$sig_x

      z_u[pp, ] <- results[[idx]]$z_u
      pi_u[pp, 1:K_epsilon] <- results[[idx]]$pi_u
      params_u[pp, 1:K_epsilon, ] <- results[[idx]]$params_u

      x_nz[pp, ] <- results[[idx]]$x_nz
    }
  }else if (!parallel){
    for (pp in (q + 1):(q + p)) {
      if (p == 0) {
        next
      }
      obj_tmp <- obj[, c(1, (1 + pp))]
      class(obj_tmp)[[1]] <- 'uninc'
      control_uni <- control
      control_uni$nsamp <- niter_uni
      control_uni$nthin <- 1
      control_uni$nburn <- 0
      if (progress){
        msg <- sprintf("preprocessing variable %s", names(obj_tmp)[2])
        cat(sprintf("\r%-*s", 50, msg))  # left-align + pad with spaces
        utils::flush.console()
      }
      res <- bdeconv(obj_tmp, lwr = lwr[pp], upr = upr[pp],
                     res_type = 'final', err_type = err_type, var_type = var_type,
                     control = control_uni, progress = FALSE)
      x[pp, ] <- res$x
      u[pp, ] <- res$u
      knots[pp, ] <- res$knots
      delta[pp] <- res$knots[2] - res$knots[1]
      theta_v[pp, ] <- res$theta_v

      z_x[pp - q, ] <- res$z_x
      pi_x[pp - q, 1:K_x] <- res$pi_x
      params_x[pp - q, 1:K_x, 1] <- res$mu_x
      params_x[pp - q, 1:K_x, 2] <- res$sig_x

      z_u[pp, ] <- res$z_u
      pi_u[pp, 1:K_epsilon] <- res$pi_u
      params_u[pp, 1:K_epsilon, ] <- cbind(res$p_u, res$mu_u, res$sig1_u,
                                         res$sig2_u)

      x_nz[pp, ] <- x[pp, ]
    }
  }
  rng_start_xs = apply(t(apply(x, 1, range)), 1, diff)
  for (dd in 1:d) {
    x_nz[dd, x_nz[dd, ] > max(knots[dd, ])] <- max(knots[dd, ]) -
      rng_start_xs[dd] / 1000
    x_nz[dd, x_nz[dd, ] < min(knots[dd, ])] <- min(knots[dd, ]) +
      rng_start_xs[dd] / 1000
    prop_sig_theta[dd, , ] <- corpcor::make.positive.definite(
      prop.sig.thetas.fn(theta_v[dd, ], x_nz[dd, ], m,
                         u[dd, ], 0.01, K_t, P2.mat(K_t + 1), knots[dd, ], n)
    )
  }
  if (q>0){
  for (qq in 1:q){
    prop_sig_theta_x_dens[qq, , ] <- diag(rep(.001,K_t+1))
  }
  }
  if(progress){
    cat(sprintf("\r Multivariate MCMC chain started \n"))
    utils::flush.console()
  }
  res <- mcmcmvt_cpp(nsamp, nburn, nthin, niter_u,
                     w, x, u, m, p, q,
                     knots, delta, theta_v, theta_x,
                     z_u - 1, pi_u,
                     ic, inc,
                     x_e, x_nz, theta_e,
                     z_x - 1, pi_x, params_x,
                     P_t, prop_sig_theta,prop_sig_theta_x_dens,
                     lwr, upr,
                     x_grd_res, var_grd_res, e_grd_res,
                     alpha_x,
                     a_sigmasq_x, b_sigmasq_x,
                     K_x, a_vartheta, b_vartheta,
                     var_type_c, K_t, err_type_c, K_epsilon, alpha_epsilon,
                     a_epsilon, b_epsilon, sigmasq_epsilon,
                     a_xi, b_xi, sigmasq_xi,
                     a_beta, b_beta, sigmasq_beta,
                     mu0_theta_e, sigma00_theta_e,
                     sig_tune_theta_1, sig_tune_theta_2,progress)
  u_res <- res$u_res
  corr_res <- res$corr_res
  res <- append(res[-c((length(res) - 1), length(res))],
                append(u_res, corr_res))
  res$q <- q
  res$p <- p
  res$w_var <- t(w_var)
  res$knots <- knots
  res$lwr <- lwr
  res$upr <- upr
  res$control <- control
  res$res_type <- res_type
  res$err_type <- err_type
  return(bdeconvmv_result(res, res_type))
}
