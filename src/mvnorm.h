#ifndef __MVNORM__
#define __MVNORM__

#include <RcppArmadillo.h>

arma::vec dmvnorm(arma::mat const &x, arma::rowvec const &mean,
                  arma::mat const &sigma, bool const logd = false);

arma::rowvec rmvnorm(arma::vec const &mean, arma::mat const &sigma);

#endif
