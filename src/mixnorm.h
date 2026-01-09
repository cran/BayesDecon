#ifndef __MIXNORM__
#define __MIXNORM__

#include <RcppArmadillo.h>


arma::vec fu_mixnorm(arma::vec x, double mean, arma::vec sd, arma::vec pi,
                     arma::mat params);

arma::vec dist_mixnorm(arma::vec x, double mean, arma::vec sd, arma::vec pi,
                       arma::mat params);

arma::vec dist_restricted_mixnorm(arma::vec x, arma::vec mean, arma::vec sd,
                                  arma::mat params);

arma::vec d_restricted_mixnorm(arma::vec x, arma::vec mean, arma::vec sd,
                               arma::mat params);

arma::vec d_restricted_mixnorm_sclr(double x, double mean, double sd,
                                    arma::mat params);

arma::vec d_scaled_restricted_mixnorm(arma::vec x, double mean, double sd,
                                      arma::vec pi, arma::mat params);

arma::vec r_tnorm_proposal_params_restricted_mix_norm(arma::vec params);

double diff_log_curr_to_prop_tnorm_proposal_params_restricted_mix_norm(arma::vec prop_params, arma::vec curr_params);

arma::vec r_proposal_params_restricted_mix_norm(double p_a, double p_b,
                                                double sigmasq_mu_curl,
                                                double s11, double s12,
                                                double s21, double s22);

#endif
