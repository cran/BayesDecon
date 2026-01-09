#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>


arma::vec sample(const arma::vec& x, const int& size, const bool& replace,
                 const arma::vec& probs){
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}

double sample(const arma::vec& x, const arma::vec& probs){
  arma::vec val = Rcpp::RcppArmadillo::sample(x, 1, true, probs);
  return val(0);
}

arma::vec sample(const arma::vec& x, const int& size, const bool& replace){
  arma::vec probs(x.n_elem, arma::fill::value(1.0 / x.n_elem));
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}

double sample(const arma::vec& x){
  arma::vec probs(x.n_elem, arma::fill::value(1.0 / x.n_elem));
  arma::vec val = Rcpp::RcppArmadillo::sample(x, 1, true, probs);
  return val(0);
}

arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace,
                  const arma::vec& probs){
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}

arma::uword sample(const arma::uvec& x, const arma::vec& probs){
  arma::uvec val = Rcpp::RcppArmadillo::sample(x, 1, true, probs);
  return val(0);
}

arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace){
  arma::vec probs(x.n_elem, arma::fill::value(1.0 / x.n_elem));
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}

arma::uword sample(const arma::uvec& x){
  arma::vec probs(x.n_elem, arma::fill::value(1.0 / x.n_elem));
  arma::uvec val = Rcpp::RcppArmadillo::sample(x, 1, true, probs);
  return val(0);
}

arma::uvec sample(const arma::uword x, const int& size, const bool& replace,
                  const arma::vec& probs) {
  arma::uvec x_lin = arma::linspace<arma::uvec>(1, x, x);
  arma::vec probs1(probs.n_elem);
  // If all probabilities are zero, use uniform distribution
  if (arma::all(probs == 0)) {
    probs1.fill(1.0 / probs.n_elem);
  } else{
    probs1 = probs;
  }
  return sample(x_lin, size, replace, probs1);
}
