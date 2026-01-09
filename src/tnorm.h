#ifndef __TNORM__
#define __TNORM__

#include <math.h>
#include <RcppArmadillo.h>


double dtnorm_sclr(double x, double mean, double sd, double lwr, double upr,
                   bool log);

arma::vec dtnorm(arma::vec x, double mean, double sd, double lwr, double upr,
                 bool log);

arma::vec dtnorm(arma::vec x, arma::vec mean, double sd, double lwr, double upr,
                 bool log);
arma::vec dtnorm(arma::vec x, arma::vec mean, arma::vec sd, double lwr, double upr,
                 bool log);

arma::vec ptnorm(arma::vec q, double mean, double sd, double lwr, double upr,
                 bool lower_tail, bool log_p);
arma::vec ptnorm(arma::vec q, arma::vec mean, arma::vec sd, double lwr, double upr,
                 bool lower_tail, bool log_p);

arma::vec rtnorm(arma::uword n, double mean, double sd, double lwr,
                 double upr);

double rtnorm(double mean, double sd, double lwr, double upr);

arma::vec rtnorm(arma::vec mean, double sd, double lwr, double upr);

arma::vec dist_mixtnorm(arma::vec x, arma::vec pi, arma::vec mean, arma::vec sd,
                        double lwr, double upr);

arma::vec dens_mixtnorm(arma::vec x, arma::vec pi, arma::vec mean, arma::vec sd,
                        double lwr, double upr);

arma::vec rmvtnorm(arma::vec alpha, arma::vec beta, arma::mat cov);

arma::vec rmvtnorm_mono(arma::vec mean, arma::mat cov);

#endif
