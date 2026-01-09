#ifndef __SAMPLE__
#define __SAMPLE__

#include <RcppArmadillo.h>


arma::vec sample(const arma::vec& x, const int& size, const bool& replace,
                 const arma::vec& probs);

double sample(const arma::vec& x, const arma::vec& probs);

arma::vec sample(const arma::vec& x, const int& size, const bool& replace);

double sample(const arma::vec& x);

arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace,
                  const arma::vec& probs);

arma::uword sample(const arma::uvec& x, const arma::vec& probs);

arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace);

arma::uword sample(const arma::uvec& x);

arma::uvec sample(const arma::uword x, const int& size, const bool& replace,
                  const arma::vec& probs);

#endif
