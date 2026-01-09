#include "mvnorm.h"

// From https://gallery.rcpp.org/articles/dmvnorm_arma/
// By Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnorm(arma::mat const &x, arma::rowvec const &mean,
                   arma::mat const &sigma, bool const logd) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

arma::rowvec rmvnorm(arma::vec const &mean, arma::mat const &sigma) {
  arma::uword n = 1;
  arma::uword d = sigma.n_cols;
  arma::mat q = arma::randn(n, d);
  return arma::repmat(mean, 1, n).t() + q * arma::chol(sigma);
}
