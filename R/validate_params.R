validate_sclr_params <- function(nsamp, nburn, nthin, niter_u,
                                 x_grd_res, var_grd_res, e_grd_res) {
  if (nsamp < 1) {
    stop('`nsamp` must be greater than 0.')
  }
  if (nburn < 0) {
    stop('`nburn` must be nonnegative.')
  }
  if (nthin < 1) {
    stop('`nthin` must be greater than 0.')
  }
  if (niter_u < 1) {
    stop('`nter_u` must be greater than 0.')
  }
  if (x_grd_res < 1) {
    stop('`x_grd_res` must be greater than 0.')
  }
  if (var_grd_res < 1) {
    stop('`var_grd_res` must be greater than 0.')
  }
  if (e_grd_res < 1) {
    stop('`e_grd_res` must be greater than 0.')
  }
}
