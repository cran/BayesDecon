#ifndef __GETCLOSE__
#define __GETCLOSE__

#include <RcppArmadillo.h>

arma::uvec get_close(arma::rowvec a, arma::rowvec b) {
    arma::uvec out(a.n_elem);
    arma::mat  tmp(a.n_elem, b.n_elem);
    for (arma::uword aa = 0; aa < a.n_elem; aa++) {
        for (arma::uword bb = 0; bb < b.n_elem; bb++) {
        tmp(aa, bb) = abs(a(aa) - b(bb));
        }
        out(aa) = tmp.row(aa).index_min();
    }
    return out;
}

#endif
