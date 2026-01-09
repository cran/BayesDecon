// #include "formcorr.h"
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat formcorr(arma::vec b, arma::vec theta) {
  arma::uword p = b.n_elem + 1;
  arma::mat V(p, p, arma::fill::zeros);
  V(0, 0) = 1;
  V(1, 0) = b(0);
  V(1, 1) = sqrt(1.0 - b(0) * b(0));
  if (theta.n_elem > 0) {
    for (arma::uword pp = 3; pp <= p; pp++) {
      arma::uword q1pp = 1 + (pp - 3) * (pp - 2) / 2;
      arma::uword q2pp = (pp - 2) * (pp - 1) / 2;
      V(pp - 1, 0) = b(pp - 2) * sin(theta(q1pp - 1));
      if (pp > 3) {
        for (arma::uword kk = 2; kk <= (pp - 2); kk++) {
          arma::vec cos_(q1pp + kk - 2 - q1pp + 1, arma::fill::zeros);
          for (arma::uword jj = q1pp; jj <= q1pp + kk - 2; jj++) {
            cos_(jj - q1pp) = cos(theta(jj - 1));
          }
          V(pp - 1, kk - 1) = b(pp - 2) * arma::prod(cos_) *
            sin(theta(q1pp + kk - 2));
        }
      }
      arma::vec cos_(q2pp - q1pp + 1, arma::fill::zeros);
      for (arma::uword jj = q1pp; jj <= q2pp; jj++) {
        cos_(jj - q1pp) = cos(theta(jj - 1));
      }
      V(pp - 1, pp - 2) = b(pp - 2) * arma::prod(cos_);
      V(pp - 1, pp - 1) = sqrt(1.0 - b(pp - 2) * b(pp - 2));
    }
  }
  return V * V.t();
}
