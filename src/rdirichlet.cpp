#include "rdirichlet.h"

arma::vec rdirichlet(arma::vec alpha) {
  arma::vec y(alpha.n_elem);
  for (arma::uword i = 0; i < alpha.n_elem; i++) {
    y(i) = R::rgamma(alpha(i), 1);
  }
  return y / sum(y);
}
