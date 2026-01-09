#include <math.h>
#include <RcppArmadillo.h>
#include "B_basis.h"
#include "tapply.h"
#include "mvnorm.h"
#include "tnorm.h"
#include "sample.h"
#include "rdirichlet.h"
#include "mixnorm.h"
#include "var_e_fn.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List mcmcunireg_cpp(arma::uword nsamp, arma::uword nburn,
                          arma::uword nthin, arma::uword niter_u,
                          arma::vec ws, arma::vec xs, arma::vec us,
                          arma::uvec mis, arma::uvec inds,
                          double xs_lwr, double xs_upr, arma::vec srvywgt,
                          bool type_v, arma::uword K_t,
                          arma::vec knots_t, arma::vec thetas, arma::mat P_t,
                          arma::uword x_grd_res, arma::uword e_grd_res,
                          arma::uword var_grd_res,
                          arma::uword K_x,
                          arma::vec pi_xs, arma::uvec z_xs,
                          arma::vec mu_xs, arma::vec sigmasq_xs,
                          arma::uword type_u, arma::uword K_epsilon,
                          double alpha_xs, double alpha_epsilon,
                          double mu0_xs, double sigmasq0_xs,
                          double a_sigmasq_xs, double b_sigmasq_xs,
                          double a_epsilon, double b_epsilon, double sigmasq_epsilon,
                          double sigmasq_vartheta, double a_vartheta, double b_vartheta,
                          arma::mat prop_sig_thetas, double sig_tune_thetas_1,
                          double sig_tune_thetas_2,
                          bool use_wgt, bool show_progress) {
  // ~~~~~~~~~~~~~ //
  // Prepare input //
  // ~~~~~~~~~~~~~ //
  arma::uword n = mis.n_elem;
  arma::uword N = sum(mis);
  arma::uvec inds1(n);
  arma::uvec inds2(n);
  inds1(0) = 1;
  inds2(0) = mis(0);
  for (arma::uword nn = 1; nn < n; nn++) {
    inds1(nn) = inds1(nn - 1) + mis(nn - 1);
    inds2(nn) = inds1(nn) + mis(nn) - 1;
  }
  inds1 -= 1;
  inds2 -= 1;
  arma::uvec ic = arma::find(ws != 0);
  arma::uvec inc = arma::find(ws == 0);

  arma::vec xs_grid = arma::linspace<arma::vec>(xs_lwr, xs_upr, x_grd_res);
  arma::vec es_grid = arma::linspace<arma::vec>(-3.0, 3.0, e_grd_res); // not needed
  arma::vec var_grid = arma::linspace<arma::vec>(xs_lwr, xs_upr, var_grd_res); // not needed

  arma::mat B_basis_store = B_basis(xs_grid, knots_t);

  arma::vec current_thetas = thetas;
  arma::vec vars = B_basis(xs, knots_t) * exp(current_thetas);
  arma::vec current_vars = vars;
  arma::vec proposed_vars = vars;

  // ~~~~~~~~~~~~~~~~ //
  // For progress bar //
  // ~~~~~~~~~~~~~~~~ //
  arma::uword next_update = 5;

  // ~~~~~~~~~~~~~~~~~~~~~~~ //
  // Storage for MCMC output //
  // ~~~~~~~~~~~~~~~~~~~~~~~ //
  arma::mat d_ordinates_xs(n, K_x);

  arma::vec proposed_xs = xs;
  arma::vec current_xs = xs;
  arma::vec proposed_us = us;
  arma::vec current_us = us;

  arma::vec current_likelihood(n, arma::fill::ones);
  arma::vec proposed_likelihood(n, arma::fill::ones);

  arma::uvec z_us(N, arma::fill::ones);
  arma::vec pi_us(K_epsilon, arma::fill::value(1.0 / K_epsilon));
  arma::mat params_us(K_epsilon, 4);
  arma::rowvec params_us_ = {0.5, 0, 1, 1};
  for (arma::uword kk = 0; kk < K_epsilon; kk++) {
    params_us.row(kk) = params_us_;
  }
  double var_es = 1;

  arma::mat thetas_mcmc(nsamp, thetas.n_elem, arma::fill::zeros);
  arma::vec var_es_mcmc(nsamp);
  arma::mat us_mcmc(nsamp, N);
  arma::vec k_us_mcmc(nsamp);
  arma::umat z_us_mcmc(nsamp, N);
  arma::mat pi_us_mcmc(nsamp, K_epsilon);
  arma::mat p_us_mcmc(nsamp, K_epsilon);
  arma::mat mu_us_mcmc(nsamp, K_epsilon);
  arma::mat sig1_us_mcmc(nsamp, K_epsilon);
  arma::mat sig2_us_mcmc(nsamp, K_epsilon);
  arma::mat xs_mcmc(nsamp, n);
  arma::umat z_xs_mcmc(nsamp, n);
  arma::mat pi_xs_mcmc(nsamp, K_x);
  arma::mat mu_xs_mcmc(nsamp, K_x);
  arma::mat sig_xs_mcmc(nsamp, K_x);
  // ~~~~~~~~~~~~ //
  // Conduct MCMC //
  // ~~~~~~~~~~~~ //
  arma::uword niter = nsamp * nthin + nburn - nthin + 1;
  for (arma::uword it = 0; it < niter; it++) {
    // Update z_xs
    for (arma::uword kk = 0; kk < K_x; kk++) {
      d_ordinates_xs.col(kk) = dtnorm(xs, mu_xs(kk), sqrt(sigmasq_xs(kk)),
                                      xs_lwr, xs_upr, false);
    }
    d_ordinates_xs.replace(arma::datum::nan, 0);
    d_ordinates_xs.elem(arma::find_nonfinite(d_ordinates_xs)).zeros();
    for (arma::uword ii = 0; ii < n; ii++) {
      arma::uvec samp_res = sample(K_x, 1, false,
                                   pi_xs % d_ordinates_xs.row(ii).t());
      z_xs(ii) = samp_res(0);
    }
    // Update cluster probabilities
    arma::uvec n_kk_xs = arma::hist(z_xs,
                                    arma::linspace<arma::uvec>(1, K_x,
                                                               K_x));
    pi_xs = rdirichlet(arma::conv_to<arma::vec>::from(n_kk_xs) + alpha_xs /
      K_x);
    if(((it >= nburn/2) & (it < nburn)) & (K_x > 5)){
      arma::uvec zero_el_cluster = arma::find(n_kk_xs == 0);
      if (!(zero_el_cluster.n_elem == 0)){
        arma::uvec ind_to_remove({});
        if (K_x - zero_el_cluster.n_elem >5){
          ind_to_remove = zero_el_cluster;
        } else {
          ind_to_remove = zero_el_cluster(arma::span(0,K_x-6));
        }
        K_x = K_x - ind_to_remove.n_elem;
        for(arma::uword elem=0; elem < ind_to_remove.n_elem-1; elem++){
          arma::uvec temp_vec = arma::find((z_xs - 1 >= ind_to_remove(elem)) && (z_xs - 1 < ind_to_remove(elem + 1)));
          z_xs.elem(temp_vec) -= (elem+1);
        }
        arma::uvec temp_vec = arma::find(z_xs - 1 >= ind_to_remove(ind_to_remove.n_elem-1));
        z_xs.elem(temp_vec) -= (ind_to_remove.n_elem);

        pi_xs.shed_rows(ind_to_remove);
        mu_xs.shed_rows(ind_to_remove);
        sigmasq_xs.shed_rows(ind_to_remove);
        d_ordinates_xs.resize(n, K_x);
        pi_xs_mcmc.resize(nsamp, K_x);
        mu_xs_mcmc.resize(nsamp, K_x);
        sig_xs_mcmc.resize(nsamp, K_x);
      }
    }
    // Update mu_xs and sigmasq_xs
    arma::vec xs_trans(n);
    for (arma::uword nn = 0; nn < n; nn++) {
      double mu_ = mu_xs(z_xs(nn) - 1);
      double sig_ = sqrt(sigmasq_xs(z_xs(nn) - 1));
      xs_trans(nn) = mu_ + sig_ * R::qnorm(
        (R::pnorm((xs(nn) - mu_) / sig_, 0, 1, true, false) -
         R::pnorm((xs_lwr - mu_) / sig_, 0, 1, true, false)) /
        (R::pnorm((xs_upr - mu_) / sig_, 0, 1, true, false) -
         R::pnorm((xs_lwr - mu_) / sig_, 0, 1, true, false)),
        0, 1, true, false);
    }
    xs_trans.clamp(xs_lwr - 10, xs_upr + 10);

    for (arma::uword kk = 0; kk < K_x; kk++) {
      arma::uvec tmp = arma::find((z_xs - 1) == kk);
      arma::vec xspool = xs_trans(tmp);
      double sigmasq_tmp = 1.0 / (n_kk_xs(kk) / sigmasq_xs(kk) +
        1.0 / sigmasq0_xs);
      double mu_tmp = (sum(xspool) / sigmasq_xs(kk) + mu0_xs / sigmasq0_xs) *
        sigmasq_tmp;

      mu_xs(kk) = R::rnorm(mu_tmp, sqrt(sigmasq_tmp));
      double post_a_sigmasq_xs = a_sigmasq_xs + xspool.n_elem / 2.0;
      double post_b_sigmasq_xs = b_sigmasq_xs + sum(pow(xspool - mu_xs(kk), 2))
        / 2.0;
      sigmasq_xs(kk) = 1.0 / R::rgamma(post_a_sigmasq_xs,
                 1.0 / post_b_sigmasq_xs);
    }
    // Update xs and us
    arma::vec proposed_xs(n);
    for (arma::uword nn = 0; nn < n; nn++) {
      arma::vec val = rtnorm(1, current_xs(nn), 0.1, xs_lwr, xs_upr);
      proposed_xs(nn) = val(0);
    }
    arma::mat tempmat(n, x_grd_res);
    for (arma::uword rr = 0; rr < n; rr++) {
      for (arma::uword cc = 0; cc < x_grd_res; cc++) {
        tempmat(rr, cc) = abs(proposed_xs(rr) - xs_grid(cc));
      }
    }
    arma::uvec close_ind(n);
    for (arma::uword rr = 0; rr < n; rr++) {
      close_ind(rr) = tempmat.row(rr).index_min();
    }
    proposed_vars = B_basis_store.rows(close_ind) * exp(thetas);

    arma::vec proposed_prior(n);
    arma::vec current_prior(n);
    arma::vec current_sds_rep(N);
    arma::vec proposed_sds_rep(N);
    arma::uword us_idx = 0;
    for (arma::uword nn = 0; nn < n; nn++) {
      double mu_ = mu_xs(z_xs(nn) - 1);
      double sig_ = sqrt(sigmasq_xs(z_xs(nn) - 1));
      proposed_prior(nn) = dtnorm_sclr(proposed_xs(nn), mu_, sig_, xs_lwr,
                                       xs_upr, false);
      current_prior(nn) = dtnorm_sclr(current_xs(nn), mu_, sig_, xs_lwr, xs_upr,
                                      false);
    }
    proposed_us = ws - proposed_xs(inds);
    current_sds_rep = sqrt(current_vars(inds));
    proposed_sds_rep = sqrt(proposed_vars(inds));

    arma::uword k_us;
    arma::vec temp_current_us_likelihood(N);
    arma::vec temp_proposed_us_likelihood(N);
    if (type_u == 0) {
      k_us = max(z_us);
      temp_current_us_likelihood = fu_mixnorm(current_us, 0,
                                              current_sds_rep,
                                              pi_us(arma::span(0,
                                                              k_us - 1
                                                              )),
                                              params_us.rows(
                                                arma::span(0, k_us - 1)
                                                ));
      temp_proposed_us_likelihood = fu_mixnorm(proposed_us, 0,
                                              proposed_sds_rep,
                                              pi_us(arma::span(0,
                                                                k_us - 1
                                                                )),
                                              params_us.rows(
                                                arma::span(0, k_us - 1)
                                                ));
    } else if (type_u == 1) {
      for (arma::uword NN = 0; NN < N; NN++) {
        temp_current_us_likelihood(NN) = R::dnorm(current_us(NN), 0, current_sds_rep(NN), false);
        temp_proposed_us_likelihood(NN) = R::dnorm(proposed_us(NN), 0, proposed_sds_rep(NN), false);
      }
    } else if (type_u == 2) {
      for (arma::uword NN = 0; NN < N; NN++) {
        temp_current_us_likelihood(NN)  = 0.5 * R::dexp(abs(current_us(NN)), sqrt(2.0) / current_sds_rep(NN), false);
        temp_proposed_us_likelihood(NN) = 0.5 * R::dexp(abs(proposed_us(NN)), sqrt(2.0) / proposed_sds_rep(NN), false);
      }
    }
    current_likelihood = tapply_prod(temp_current_us_likelihood, mis);
    proposed_likelihood = tapply_prod(temp_proposed_us_likelihood, mis);

    arma::vec mh_ratio = proposed_prior % proposed_likelihood;
    mh_ratio /= current_prior % current_likelihood;
    for (arma::uword nn = 0; nn < n; nn++) {
      mh_ratio(nn) *= dtnorm_sclr(current_xs(nn), proposed_xs(nn), 0.1,
                                  xs_lwr, xs_upr, false);
      mh_ratio(nn) /= dtnorm_sclr(proposed_xs(nn), current_xs(nn), 0.1,
                                  xs_lwr, xs_upr, false);
    }
    mh_ratio.replace(arma::datum::nan, 0);
    arma::vec u(n, arma::fill::randu);
    arma::uvec inds_to_replace = arma::find(u < mh_ratio);
    current_xs(inds_to_replace) = proposed_xs(inds_to_replace);
    xs(inds_to_replace) = current_xs(inds_to_replace);
    current_vars(inds_to_replace) = proposed_vars(inds_to_replace);
    vars(inds_to_replace) = current_vars(inds_to_replace);
    current_us = ws - xs(inds);
    us = current_us;
    // Update thetas
    if (type_v) {
      arma::mat theta_cov(K_t + 1, K_t + 1);
      theta_cov.diag() = arma::vec(K_t + 1, arma::fill::value(sig_tune_thetas_1));
      theta_cov += sig_tune_thetas_2 * prop_sig_thetas;
      arma::vec proposed_thetas = rmvnorm(current_thetas, theta_cov).t();
      for (arma::uword rr = 0; rr < n; rr++) {
        for (arma::uword cc = 0; cc < x_grd_res; cc++) {
          tempmat(rr, cc) = abs(xs(rr) - xs_grid(cc));
        }
      }
      for (arma::uword rr = 0; rr < n; rr++) {
        close_ind(rr) = tempmat.row(rr).index_min();
      }
      proposed_vars = B_basis_store.rows(close_ind) * exp(proposed_thetas);

      arma::vec current_log_prior = -current_thetas.t() * P_t * current_thetas;
      current_log_prior /= (2.0 * sigmasq_vartheta);
      arma::vec proposed_log_prior = -proposed_thetas.t() * P_t * proposed_thetas;
      proposed_log_prior /= (2.0 * sigmasq_vartheta);
      arma::vec mean_(N, arma::fill::zeros);
      us_idx = 0;
      for (arma::uword nn = 0; nn < n; nn++) {
        for (arma::uword mm = 0; mm < mis(nn); mm++) {
          current_sds_rep(us_idx) = sqrt(current_vars(nn));
          proposed_sds_rep(us_idx) = sqrt(proposed_vars(nn));
          us_idx++;
        }
      }
      arma::vec temp_current_likelihood(N);
      arma::vec temp_proposed_likelihood(N);
      if (type_u == 0) {
        temp_current_likelihood = d_restricted_mixnorm(us, mean_,
                                                       current_sds_rep,
                                                       params_us.rows(
                                                         z_us - 1
                                                       ));
        temp_proposed_likelihood = d_restricted_mixnorm(us, mean_,
                                                        proposed_sds_rep,
                                                        params_us.rows(
                                                          z_us - 1
                                                        ));
      } else if (type_u == 1) {
        for (arma::uword NN = 0; NN < N; NN++) {
          temp_current_likelihood(NN) = R::dnorm(us(NN), 0, current_sds_rep(NN), false);
          temp_proposed_likelihood(NN) = R::dnorm(us(NN), 0, proposed_sds_rep(NN), false);
        }
      } else if (type_u == 2) {
        for (arma::uword NN = 0; NN < N; NN++) {
          temp_current_likelihood(NN)  = 0.5 * R::dexp(abs(us(NN)), sqrt(2.0) / current_sds_rep(NN), false);
          temp_proposed_likelihood(NN) = 0.5 * R::dexp(abs(us(NN)), sqrt(2.0) / proposed_sds_rep(NN), false);
        }
      }
      double current_log_likelihood = sum(log(temp_current_likelihood));
      double proposed_log_likelihood = sum(log(temp_proposed_likelihood));
      double log_mh_ratio = proposed_log_prior(0) + proposed_log_likelihood -
        current_log_likelihood - current_log_prior(0);
      if (log_mh_ratio == arma::datum::nan) {
        log_mh_ratio = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_ratio) {
        current_thetas = proposed_thetas;
        thetas = current_thetas;
        current_vars = proposed_vars;
        vars = current_vars;
      }
    }
    // Update sigmasq_vartheta
    arma::vec sigmasq_vartheta_rate = b_vartheta + thetas.t() * P_t * thetas / 2.0; // changed /2.0
    sigmasq_vartheta = 1.0 / R::rgamma(a_vartheta + (K_t + 1) / 2.0, 1.0 / sigmasq_vartheta_rate(0));

    if (type_u == 0) {
      // Update z_us
      for (arma::uword NN = 0; NN < N; NN++) {
        arma::vec prob_us = pi_us %
          d_restricted_mixnorm(arma::vec(K_epsilon, arma::fill::value(us(NN))),
                              arma::vec(K_epsilon, arma::fill::zeros),
                              arma::vec(K_epsilon,
                                        arma::fill::value(sqrt(vars(inds(NN))))),
                              params_us);
        if (sum(prob_us) == 0) {
          prob_us = arma::vec(K_epsilon, arma::fill::value(1.0 / K_epsilon));
        }
        arma::uvec samp_ = sample(K_epsilon, 1, false, prob_us); // changed false from true
        z_us(NN) = samp_(0);
      }
      // Update cluster probabilites
      arma::uvec n_kk_us = arma::hist(z_us,
                                      arma::linspace<arma::uvec>(1, K_epsilon,
                                                                 K_epsilon));
      pi_us = rdirichlet(arma::conv_to<arma::vec>::from(n_kk_us) + alpha_epsilon /
        K_epsilon);
      if(((it >= nburn/2) & (it < nburn)) & (K_epsilon > 5)){
        arma::uvec zero_el_cluster = arma::find(n_kk_us == 0);
        if (!(zero_el_cluster.n_elem == 0)){
          arma::uvec ind_to_remove({});
          if (K_epsilon - zero_el_cluster.n_elem >5){
            ind_to_remove = zero_el_cluster;
          } else {
            ind_to_remove = zero_el_cluster(arma::span(0,K_epsilon-6));
          }
          K_epsilon = K_epsilon - ind_to_remove.n_elem;
          for(arma::uword elem=0; elem < ind_to_remove.n_elem-1; elem++){
            arma::uvec temp_vec = arma::find((z_us - 1 >= ind_to_remove(elem)) && (z_us - 1 < ind_to_remove(elem + 1)));
            z_us.elem(temp_vec) -= (elem+1);
          }
          arma::uvec temp_vec = arma::find(z_us - 1 >= ind_to_remove(ind_to_remove.n_elem-1));
          z_us.elem(temp_vec) -= (ind_to_remove.n_elem);

          params_us.shed_rows(ind_to_remove);
          pi_us.shed_rows(ind_to_remove);
          pi_us_mcmc.resize(nsamp, K_epsilon);
          p_us_mcmc.resize(nsamp, K_epsilon);
          mu_us_mcmc.resize(nsamp, K_epsilon);
          sig1_us_mcmc.resize(nsamp, K_epsilon);
          sig2_us_mcmc.resize(nsamp, K_epsilon);
        }
      }
      // Update params_us
      k_us = K_epsilon;
      if (it > 2000) {
        niter_u = 1;
      }
      for (arma::uword rr = 0; rr < niter_u; rr++) {
        for (arma::uword kk = 0; kk < k_us; kk++) {
          arma::uvec temp = arma::find(z_us - 1 == kk);
          arma::vec uspool = us(temp);
          arma::vec varspool = vars(inds(temp));
          arma::vec proposed_params_us =
            r_tnorm_proposal_params_restricted_mix_norm(params_us.row(kk).t());
          double proposed_log_prior = R::dnorm(proposed_params_us(1), 0,
                                               sqrt(sigmasq_epsilon), true) -
                                                 (a_epsilon + 1) * (log(proposed_params_us(2)) +
                                                 log(proposed_params_us(3))) - b_epsilon * (1.0 / proposed_params_us(2) +
                                                 1.0 / proposed_params_us(3));
          double current_log_prior = R::dnorm(params_us(kk, 1), 0,
                                              sqrt(sigmasq_epsilon), true) -
                                                (a_epsilon + 1) * (log(params_us(kk, 2)) + log(params_us(kk, 3))) -
                                                b_epsilon * (1.0 / params_us(kk, 2) + 1.0 / params_us(kk, 3));
          arma::vec temp_proposed_log_likelihood = log(
            d_restricted_mixnorm(uspool,
                                arma::vec(uspool.n_elem, arma::fill::zeros),
                                sqrt(varspool),
                                proposed_params_us.t())
          );
          arma::vec temp_current_log_likelihood = log(
            d_restricted_mixnorm(uspool,
                                arma::vec(uspool.n_elem, arma::fill::zeros),
                                sqrt(varspool),
                                params_us.row(kk))
          );
          temp_proposed_log_likelihood.replace(arma::datum::nan, 0);
          temp_current_log_likelihood.replace(arma::datum::nan, 0);
          double proposed_log_likelihood = sum(temp_proposed_log_likelihood);
          double current_log_likelihood = sum(temp_current_log_likelihood);

          double log_acc_prob = proposed_log_prior + proposed_log_likelihood -
            current_log_prior - current_log_likelihood;
          log_acc_prob += diff_log_curr_to_prop_tnorm_proposal_params_restricted_mix_norm(proposed_params_us,params_us.row(kk).t());//changed
          if (log(R::runif(0, 1)) < log_acc_prob) {
            params_us.row(kk) = proposed_params_us.t();
          }
        }
      }
      if (k_us < K_epsilon) {
        for (arma::uword kk = k_us; kk < K_epsilon; kk++) {
          params_us.row(kk) = r_proposal_params_restricted_mix_norm(1, 1, 3, 3, 3, 3, 3).t();
        }
      }

       var_es = var_e_fn(pi_us(arma::span(0, k_us - 1)),
                           params_us.rows(arma::span(0, k_us - 1)));
       params_us(arma::span(0, k_us - 1), 1) /= sqrt(var_es);
       params_us(arma::span(0, k_us - 1), arma::span(2, 3)) /= var_es;
    }


    if ((it >= nburn) & ((it - nburn) % nthin == 0)) {
      arma::uword idx = (it - nburn) / nthin;
      thetas_mcmc.row(idx) = thetas.t();
      xs_mcmc.row(idx) = xs.t();
      z_xs_mcmc.row(idx) = z_xs.t();
      pi_xs_mcmc.row(idx) = pi_xs.t();
      mu_xs_mcmc.row(idx) = mu_xs.t();
      sig_xs_mcmc.row(idx) = sqrt(sigmasq_xs.t());
      us_mcmc.row(idx) = us.t();
      if (type_u == 0) {
        var_es = var_e_fn(pi_us(arma::span(0, k_us - 1)),
                          params_us.rows(arma::span(0, k_us - 1)));
        k_us_mcmc(idx) = max(z_us);
        z_us_mcmc.row(idx) = z_us.t();
        pi_us_mcmc.row(idx) = pi_us.t();
        p_us_mcmc.row(idx) = params_us.col(0).t();
        mu_us_mcmc.row(idx) = params_us.col(1).t();
        sig1_us_mcmc.row(idx) = sqrt(params_us.col(2).t());
        sig2_us_mcmc.row(idx) = sqrt(params_us.col(3).t());
      } else {
        var_es = 1;
      }
      if (!type_v) {
        var_es *= 1.0 / R::rgamma(1 + 0.5 * N, 1.0 / (1 + 0.5 * arma::accu(us % us)));
      }
      var_es_mcmc(idx) = var_es;
    }
    R_CheckUserInterrupt();
    // ~~~~~~~~~~~~~~~~ //
    // For progress bar //
    // ~~~~~~~~~~~~~~~~ //

    if (show_progress){
      arma::uword percent = (100 * it) / (niter-1);
      if (percent >= next_update || it == (niter-1)) {
        arma::uword width = 50;
        arma::uword pos = (width * it) / (niter-1);

        std::string bar = "[";
        for (arma::uword jj = 0; jj < width; ++jj) {
          if (jj < pos) bar += "=";
          else if (jj == pos) bar += ">";
          else bar += " ";
        }
        bar += "]";

        Rcpp::Rcout << "\r" << bar << " " << percent << "%";
        Rcpp::Rcout.flush();

        next_update += 5;
      }
    }
  }
  return Rcpp::List::create(Rcpp::_("x") = xs_mcmc,
                            Rcpp::_("z_x") = z_xs_mcmc,
                            Rcpp::_("pi_x") = pi_xs_mcmc,
                            Rcpp::_("mu_x") = mu_xs_mcmc,
                            Rcpp::_("sig_x") = sig_xs_mcmc,
                            Rcpp::_("u") = us_mcmc,
                            Rcpp::_("k_u") = k_us_mcmc,
                            Rcpp::_("z_u") = z_us_mcmc,
                            Rcpp::_("pi_u") = pi_us_mcmc,
                            Rcpp::_("p_u") = p_us_mcmc,
                            Rcpp::_("mu_u") = mu_us_mcmc,
                            Rcpp::_("sig1_u") = sig1_us_mcmc,
                            Rcpp::_("sig2_u") = sig2_us_mcmc,
                            Rcpp::_("theta_v") = thetas_mcmc,
                            Rcpp::_("knots") = knots_t,
                            Rcpp::_("var_e") = var_es_mcmc);
}
