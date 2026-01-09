#include "B_basis.h"

arma::mat B_basis(arma::vec x, arma::vec knots) {
  arma::uword n = x.n_elem;
  arma::uword K = knots.n_elem;
  arma::mat B(n, K + 1, arma::fill::zeros);
  for (arma::uword j = 0; j < (K - 1); j++) {
    arma::uvec idxs = arma::find((x >= knots(j)) && (x <= knots(j + 1)));
    arma::vec act_x = x(idxs);
    arma::vec rsc_x = (act_x - knots(j)) / (knots(j + 1) - knots(j));
    arma::uvec j_(1, arma::fill::value(j));
    B.submat(idxs, j_) = 0.5 * pow(1.0 - rsc_x, 2);
    B.submat(idxs, j_ + 1) = -pow(rsc_x, 2) + rsc_x + 0.5;
    B.submat(idxs, j_ + 2) = 0.5 * pow(rsc_x, 2);
  }
  return B;
}

arma::vec B_basis_density_coeffs(arma::vec thetadens, double delta,
                                 arma::uword K) {
  arma::vec coefs = exp(thetadens) /
    ((delta / 6) *
      (exp(thetadens(0)) +
      5 * exp(thetadens(1)) +
      6 * sum(exp(thetadens.rows(2, K - 2))) +
      5 * exp(thetadens(K - 1)) +
      exp(thetadens(K))));
  return coefs;
}

arma::mat B_basis_cum_int(arma::vec x, arma::vec knots) {
  double delta = knots(1) - knots(0);
  arma::uword n = x.n_elem;
  arma::uword K = knots.n_elem;
  arma::mat B(n, K + 1, arma::fill::zeros);
  for (arma::uword j = 0; j < (K - 1); j++) {
    arma::uvec idxs = arma::find(x >= knots(j));
    arma::vec act_x = x(idxs);
    arma::vec rsc_x = (act_x - knots(j)) / (knots(j + 1) - knots(j));
    rsc_x.clamp(-arma::datum::inf, 1);
    arma::uvec j_(1, arma::fill::value(j));

    B.submat(idxs, j_) = 5 * delta / 6 + (delta / 2) *
      (rsc_x - pow(rsc_x, 2) + pow(rsc_x, 3) / 3);
    B.submat(idxs, j_ + 1) = delta / 6 + delta *
      (-pow(rsc_x, 3) / 3 + pow(rsc_x, 2) / 2 + rsc_x / 2);
    B.submat(idxs, j_ + 2) = delta * pow(rsc_x, 3) / 6;
  }
  arma::uvec idxs = arma::find(x >= knots(0));
  arma::uvec one_(1, arma::fill::ones);
  B.submat(idxs, one_ - 1) -= 5 * delta / 6;
  B.submat(idxs, one_) -= delta / 6;

  return B;
}

arma::vec dist_norm_basis(arma::vec x, arma::vec knots, arma::vec theta) {
  double delta = knots(1) - knots(0);
  arma::vec dist = B_basis_cum_int(x, knots) *
    B_basis_density_coeffs(theta, delta, knots.n_elem);
  return dist;
}

