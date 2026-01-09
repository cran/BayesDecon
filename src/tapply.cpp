#include "tapply.h"

arma::vec tapply_sum(arma::vec vals, arma::uvec mis) {
  arma::uword n = mis.n_elem;
  arma::vec out(n);
  arma::uword m_sum = 0;
  for (arma::uword nn = 0; nn < n; nn++) {
    arma::uword m = mis(nn);
    out(nn) = sum(vals(arma::span(m_sum, m_sum + m - 1)));
    m_sum += m;
  }
  return out;
}

arma::vec tapply_mean(arma::vec vals, arma::uvec mis) {
  arma::uword n = mis.n_elem;
  arma::vec out(n);
  arma::uword m_sum = 0;
  for (arma::uword nn = 0; nn < n; nn++) {
    arma::uword m = mis(nn);
    out(nn) = mean(vals(arma::span(m_sum, m_sum + m - 1)));
    m_sum += m;
  }
  return out;
}

arma::vec tapply_var(arma::vec vals, arma::uvec mis) {
  arma::uword n = mis.n_elem;
  arma::vec out(n);
  arma::uword m_sum = 0;
  for (arma::uword nn = 0; nn < n; nn++) {
    arma::uword m = mis(nn);
    out(nn) = var(vals(arma::span(m_sum, m_sum + m - 1)));
    m_sum += m;
  }
  return out;
}

arma::vec tapply_prod(arma::vec vals, arma::uvec mis) {
  arma::uword n = mis.n_elem;
  arma::vec out(n);
  arma::uword m_sum = 0;
  for (arma::uword nn = 0; nn < n; nn++) {
    arma::uword m = mis(nn);
    out(nn) = prod(vals(arma::span(m_sum, m_sum + m - 1)));
    m_sum += m;
  }
  return out;
}
