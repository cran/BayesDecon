#ifndef __TAPPLY__
#define __TAPPLY__

#include <RcppArmadillo.h>

arma::vec tapply_sum(arma::vec vec, arma::uvec mis);

arma::vec tapply_mean(arma::vec vec, arma::uvec mis);

arma::vec tapply_var(arma::vec vec, arma::uvec mis);

arma::vec tapply_prod(arma::vec vec, arma::uvec mis);

#endif
