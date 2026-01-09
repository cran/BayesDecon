#' Constructor for `zi_res` Objects
#'
#' @param obj A list as returned by bdeconv.unizi().
#' @param res_type The type of result returned from the Markov chain Monte Carlo
#'  sampler.
#'
#' @return A `zi_res` object.
#' @seealso [bdeconv()]
#' @keywords internal
bdeconvzi_result <- function(obj, res_type = c('final', 'mean', 'all')) {
  res_type <- match.arg(res_type)
  obj$var_e <- as.numeric(obj$var_e)
  obj$knots <- as.numeric(obj$knots)
  obj$k_u <- as.numeric(obj$k_u)
  obj <- structure(obj, class = c('bdeconv_result', 'uv_res', 'zi_res'))
  if (res_type == 'mean') {
    return(mean(obj))
  } else if (res_type == 'final') {
    obj <- obj[[length(obj)]]
  }
  return(obj)
}

#' Constructor for `nc_res` Objects
#'
#' @param obj A list as returned by bdeconv.uninc().
#' @param res_type The type of result returned from the Markov chain Monte Carlo
#'  sampler.
#'
#' @return An `nc_res` object.
#' @seealso [bdeconv()]
#' @keywords internal
bdeconvnc_result <- function(obj, res_type = c('final', 'mean', 'all')) {
  res_type <- match.arg(res_type)
  obj$w <- as.numeric(obj$w)
  obj$var_e <- as.numeric(obj$var_e)
  obj$knots <- as.numeric(obj$knots)
  obj$k_u <- as.numeric(obj$k_u)
  obj <- structure(obj, class = c('bdeconv_result', 'uv_res', 'nc_res'))
  if (res_type == 'mean') {
    obj <- mean(obj)
  } else if (res_type == 'final') {
    obj <- obj[[length(obj)]]
  }
  return(obj)
}



#' Constructor for `mv_res` Objects
#'
#' @param obj A list as returned by bdeconv.mvtzi().
#' @param res_type The type of result returned from the Markov chain Monte Carlo
#'  sampler.
#'
#' @return An `mv_res` object.
#' @seealso [bdeconv()]
#' @keywords internal
bdeconvmv_result <- function(obj, res_type = c('final', 'mean', 'all')) {
  res_type <- match.arg(res_type)
  obj <- structure(obj, class = c('bdeconv_result', 'mv_res'))
  if (res_type == 'mean') {
    return(mean(obj))
  } else if (res_type == 'final') {
    obj <- obj[[length(obj)]]
  }
  return(obj)
}

