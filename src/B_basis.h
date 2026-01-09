#ifndef __B_BASIS__
#define __B_BASIS__

#include <RcppArmadillo.h>

arma::mat B_basis(arma::vec x, arma::vec knots);

arma::vec B_basis_density_coeffs(arma::vec thetadens, double delta,
                                 arma::uword K);

arma::mat B_basis_cum_int(arma::vec x, arma::vec knots);

arma::vec dist_norm_basis(arma::vec x, arma::vec knots, arma::vec theta);

#endif
