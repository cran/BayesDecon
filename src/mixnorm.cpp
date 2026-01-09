#include <RcppArmadillo.h>
#include "tnorm.h"
#include "var_e_fn.h"

arma::vec fu_mixnorm(arma::vec x, double mean, arma::vec sd, arma::vec pi,
                     arma::mat params) {
  arma::vec p(1);
  arma::vec mu_curl(1);
  arma::vec sigmasq1(1);
  arma::vec sigmasq2(1);
  if (params.n_rows == 1) {
    p(0) = params(0, 0);
    mu_curl(0) = params(0, 1);
    sigmasq1(0) = params(0, 2);
    sigmasq2(0) = params(0, 3);
  } else {
    p = params.col(0);
    mu_curl = params.col(1);
    sigmasq1 = params.col(2);
    sigmasq2 = params.col(3);
  }

  arma::vec c0 = (p % p + (1 - p) % (1 - p));
  arma::vec c1 = (1 - p) / c0;
  arma::vec c2 = - p / c0;

  arma::vec mu1 = c1 % mu_curl;
  arma::vec mu2 = c2 % mu_curl;

  arma::mat y(p.n_elem, x.n_elem, arma::fill::zeros);
  for (arma::uword kk = 0; kk < p.n_elem; kk++) {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = (p(kk) * R::dnorm(x(nn), mean + sd(nn) * mu1(kk),
                        sd(nn) * sqrt(sigmasq1(kk)), false));
      double comp2 = (1.0 - p(kk)) * R::dnorm(x(nn), mean + sd(nn) * mu2(kk),
                                               sd(nn) * sqrt(sigmasq2(kk)),
                                               false);
      y(kk, nn) = pi(kk) * (comp1 + comp2);
    }
  }
  return sum(y, 0).t();
}

arma::vec dist_mixnorm(arma::vec x, double mean, arma::vec sd, arma::vec pi,
                       arma::mat params) {
  arma::vec p(1);
  arma::vec mu_curl(1);
  arma::vec sigmasq1(1);
  arma::vec sigmasq2(1);
  if (params.n_rows == 1) {
    p(0) = params(0, 0);
    mu_curl(0) = params(0, 1);
    sigmasq1(0) = params(0, 2);
    sigmasq2(0) = params(0, 3);
  } else {
    p = params.col(0);
    mu_curl = params.col(1);
    sigmasq1 = params.col(2);
    sigmasq2 = params.col(3);
  }

  arma::vec c0 = (p % p + (1 - p) % (1 - p));
  arma::vec c1 = (1 - p) / c0;
  arma::vec c2 = - p / c0;

  arma::vec mu1 = c1 % mu_curl;
  arma::vec mu2 = c2 % mu_curl;

  arma::mat y(p.n_elem, x.n_elem, arma::fill::zeros);
  for (arma::uword kk = 0; kk < p.n_elem; kk++) {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = (p(kk) * R::pnorm(x(nn), mean + sd(nn) * mu1(kk),
                        sd(nn) * sqrt(sigmasq1(kk)), true, false));
      double comp2 = (1.0 - p(kk)) * R::pnorm(x(nn), mean + sd(nn) * mu2(kk),
                      sd(nn) * sqrt(sigmasq2(kk)),
                      true, false);
      y(kk, nn) = pi(kk) * (comp1 + comp2);
    }
  }
  return sum(y, 0).t();
}

arma::vec dist_restricted_mixnorm(arma::vec x, arma::vec mean, arma::vec sd,
                               arma::mat params) {
  arma::vec p(1);
  arma::vec mu_curl(1);
  arma::vec sigmasq1(1);
  arma::vec sigmasq2(1);
  if (params.n_rows == 1) {
    p(0) = params(0, 0);
    mu_curl(0) = params(0, 1);
    sigmasq1(0) = params(0, 2);
    sigmasq2(0) = params(0, 3);
  } else {
    p = params.col(0);
    mu_curl = params.col(1);
    sigmasq1 = params.col(2);
    sigmasq2 = params.col(3);
  }

  arma::vec c0 = (p % p + (1 - p) % (1 - p));
  arma::vec c1 = (1 - p) / c0;
  arma::vec c2 = - p / c0;

  arma::vec mu1 = c1 % mu_curl;
  arma::vec mu2 = c2 % mu_curl;

  arma::vec y(x.n_elem, arma::fill::zeros);
  if (params.n_rows == 1) {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = p(0) * R::pnorm(x(nn), mean(nn) + sd(nn) * mu1(0),
                       sd(nn) * sqrt(sigmasq1(0)), true, false);
      double comp2 = (1.0 - p(0)) * R::pnorm(x(nn), mean(nn) + sd(nn) * mu2(0),
                      sd(nn) * sqrt(sigmasq2(0)), true,
                      false);
      y(nn) = comp1 + comp2;
    }
  } else {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = p(nn) * R::pnorm(x(nn), mean(nn) + sd(nn) * mu1(nn),
                       sd(nn) * sqrt(sigmasq1(nn)), true, false);
      double comp2 = (1.0 - p(nn)) * R::pnorm(x(nn), mean(nn) + sd(nn) * mu2(nn),
                      sd(nn) * sqrt(sigmasq2(nn)), true,
                      false);
      y(nn) = comp1 + comp2;
    }
  }
  return y;
}

arma::vec d_restricted_mixnorm(arma::vec x, arma::vec mean, arma::vec sd,
                               arma::mat params) {
  arma::vec p(1);
  arma::vec mu_curl(1);
  arma::vec sigmasq1(1);
  arma::vec sigmasq2(1);
  if (params.n_rows == 1) {
    p(0) = params(0, 0);
    mu_curl(0) = params(0, 1);
    sigmasq1(0) = params(0, 2);
    sigmasq2(0) = params(0, 3);
  } else {
    p = params.col(0);
    mu_curl = params.col(1);
    sigmasq1 = params.col(2);
    sigmasq2 = params.col(3);
  }

  arma::vec c0 = (p % p + (1 - p) % (1 - p));
  arma::vec c1 = (1 - p) / c0;
  arma::vec c2 = - p / c0;

  arma::vec mu1 = c1 % mu_curl;
  arma::vec mu2 = c2 % mu_curl;

  arma::vec y(x.n_elem, arma::fill::zeros);
  if (params.n_rows == 1) {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = p(0) * R::dnorm(x(nn), mean(nn) + sd(nn) * mu1(0),
                       sd(nn) * sqrt(sigmasq1(0)), false);
      double comp2 = (1.0 - p(0)) * R::dnorm(x(nn), mean(nn) + sd(nn) * mu2(0),
                      sd(nn) * sqrt(sigmasq2(0)),
                      false);
      y(nn) = comp1 + comp2;
    }
  } else {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = p(nn) * R::dnorm(x(nn), mean(nn) + sd(nn) * mu1(nn),
                       sd(nn) * sqrt(sigmasq1(nn)), false);
      double comp2 = (1.0 - p(nn)) * R::dnorm(x(nn), mean(nn) + sd(nn) * mu2(nn),
                      sd(nn) * sqrt(sigmasq2(nn)),
                      false);
      y(nn) = comp1 + comp2;
    }
  }
  return y;
}

arma::vec d_restricted_mixnorm_sclr(double x, double mean, double sd,
                                    arma::mat params) {
  arma::vec x_(params.n_rows, arma::fill::value(x));
  arma::vec mean_(params.n_rows, arma::fill::value(mean));
  arma::vec sd_(params.n_rows, arma::fill::value(sd));
  return(d_restricted_mixnorm(x_, mean_, sd_, params));
}

arma::vec d_scaled_restricted_mixnorm(arma::vec x, double mean, double sd,
                                      arma::vec pi, arma::mat params) {
  arma::vec p(1);
  arma::vec mu_curl(1);
  arma::vec sigmasq1(1);
  arma::vec sigmasq2(1);
  if (params.n_rows == 1) {
    p(0) = params(0, 0);
    mu_curl(0) = params(0, 1);
    sigmasq1(0) = params(0, 2);
    sigmasq2(0) = params(0, 3);
  } else {
    p = params.col(0);
    mu_curl = params.col(1);
    sigmasq1 = params.col(2);
    sigmasq2 = params.col(3);
  }

  arma::vec c0 = (p % p + (1 - p) % (1 - p));
  arma::vec c1 = (1 - p) / c0;
  arma::vec c2 = - p / c0;

  arma::vec mu1 = c1 % mu_curl;
  arma::vec mu2 = c2 % mu_curl;

  double sd_e = sqrt(var_e_fn(pi, params));
  arma::mat y(p.n_elem, x.n_elem, arma::fill::zeros);
  for (arma::uword kk = 0; kk < p.n_elem; kk++) {
    for (arma::uword nn = 0; nn < x.n_elem; nn++) {
      double comp1 = p(kk) * R::dnorm(x(nn), (mean + sd * mu1(kk)) / sd_e,
                       sd * sqrt(sigmasq1(kk)) / sd_e, false);
      double comp2 = (1.0 - p(kk)) * R::dnorm(x(nn),
                      (mean + sd * mu2(kk)) / sd_e,
                      sd * sqrt(sigmasq2(kk)) / sd_e,
                      false);
      y(kk, nn) = pi(kk) * (comp1 + comp2);
    }
  }
  return sum(y, 0).t();
}

arma::vec r_tnorm_proposal_params_restricted_mix_norm(arma::vec params) {
  double current_p = params(0);
  double current_mu_curl = params(1);
  double current_sigmasq1 = params(2);
  double current_sigmasq2 = params(3);

  arma::vec p_ = rtnorm(1, current_p, 0.01, 0, 1);
  double p = p_(0);

  arma::vec sigmasq1_ = rtnorm(1, current_sigmasq1, 0.1, 0, arma::datum::inf); // .1
  double sigmasq1 = sigmasq1_(0);
  arma::vec sigmasq2_ = rtnorm(1, current_sigmasq2, 0.1, 0, arma::datum::inf); // .1
  double sigmasq2 = sigmasq2_(0);
  double mu_curl = R::rnorm(current_mu_curl, 0.1); // .01

  arma::vec y(4);
  y(0) = p;
  y(1) = mu_curl;
  y(2) = sigmasq1;
  y(3) = sigmasq2;
  return y;
}

double diff_log_curr_to_prop_tnorm_proposal_params_restricted_mix_norm(arma::vec prop_params, arma::vec curr_params) {
  double current_p = curr_params(0);
  double current_sigmasq1 = curr_params(2);
  double current_sigmasq2 = curr_params(3);
  double proposed_p = prop_params(0);
  double proposed_sigmasq1 = prop_params(2);
  double proposed_sigmasq2 = prop_params(3);
  double num = 0.0;
  double den = 0.0;
  double diff_num = R::pnorm(1, current_p, .01, true, false) -
    R::pnorm(0, current_p, .01, true, false);
  double diff_den = R::pnorm(1, proposed_p, .01, true, false) -
    R::pnorm(0, proposed_p, .01, true, false);

  if (diff_num > 0 && R_finite(diff_num)) {
    num += std::log(diff_num);
  }

  if (diff_den > 0 && R_finite(diff_den)) {
    den += std::log(diff_den);
  }
  num += R::pnorm(0, current_sigmasq1, 0.1, false, true);
  den += R::pnorm(0, proposed_sigmasq1, 0.1, false, true);
  num += R::pnorm(0, current_sigmasq2, 0.1, false, true);
  den += R::pnorm(0, proposed_sigmasq2, 0.1, false, true);
  return num - den;
}

arma::vec r_proposal_params_restricted_mix_norm(double p_a, double p_b,
                                                double sigmasq_mu_curl,
                                                double s11, double s12,
                                                double s21, double s22) {
  double p = R::rbeta(p_a, p_b);

  double mu_curl = R::rnorm(0, sqrt(sigmasq_mu_curl));
  double sigmasq1 = 1.0 / R::rgamma(s11, 1.0 / s12);
  double sigmasq2 = 1.0 / R::rgamma(s21, 1.0 / s22);

  arma::vec y(4);
  y(0) = p;
  y(1) = mu_curl;
  y(2) = sigmasq1;
  y(3) = sigmasq2;

  return y;
}
