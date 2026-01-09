#include <RcppArmadillo.h>

double var_e_fn(arma::vec pi, arma::mat params) {
  arma::vec p = params.col(0);
  arma::vec mu_curl = params.col(1);
  arma::vec sigmasq1 = params.col(2);
  arma::vec sigmasq2 = params.col(3);

  arma::vec c0 = p % p + (1.0 - p) % (1.0 - p);
  arma::vec c1 = (1.0 - p) / c0;
  arma::vec c2 = -p / c0;

  arma::vec mu1 = c1 % mu_curl;
  arma::vec mu2 = c2 % mu_curl;

  arma::vec y = p % (mu1 % mu1 + sigmasq1) + (1.0 - p) % (mu2 % mu2 + sigmasq2);
  double y_ = sum(pi % y);
  return y_;
}
