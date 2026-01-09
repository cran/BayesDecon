#include <math.h>
#include <RcppArmadillo.h>
#include "B_basis.h"
#include "mvnorm.h"
#include "tnorm.h"
#include "mixnorm.h"
#include "sample.h"
#include "rdirichlet.h"
#include "var_e_fn.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List mcmcuniepi_cpp(arma::uword nsamp, arma::uword nburn,
                          arma::uword nthin, arma::uword niter_u,
                          arma::vec ws, arma::vec xs, arma::vec us,
                          arma::uvec mis, arma::uvec inds,
                          double xs_lwr, double xs_upr,
                          bool type_v, arma::uword K_t,
                          arma::vec knots_t, arma::vec thetas, arma::mat P_t,
                          arma::uword x_grd_res, arma::uword e_grd_res,
                          arma::uword var_grd_res,
                          arma::vec thetas_xs_density,
                          arma::uword type_u, arma::uword K_epsilon,
                          double a_vartheta, double b_vartheta, double sigmasq_vartheta,
                          double a_beta, double b_beta, double sigmasq_beta,
                          double a_xi, double b_xi,
                          double sigmasq_xi,
                          double alpha_epsilon, double a_epsilon, double b_epsilon,
                          double sigmasq_epsilon,
                          arma::vec mu0_thetas_eps, arma::mat sigma0_thetas_eps,
                          arma::mat prop_sig_thetas_xs_density, arma::mat prop_sig_thetas,
                          double sig_tune_thetas_1, double sig_tune_thetas_2,
                          bool show_progress)  {
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
  double range_start_xs = max(xs) - min(xs);

  arma::vec xs_grid = arma::linspace<arma::vec>(xs_lwr, xs_upr, x_grd_res);
  arma::vec es_grid = arma::linspace<arma::vec>(-3.0, 3.0, e_grd_res);
  arma::vec var_grid = arma::linspace<arma::vec>(xs_lwr, xs_upr, var_grd_res);

  double delta_t = knots_t(1) - knots_t(0);

  arma::vec thetas_eps = mu0_thetas_eps;
  arma::vec proposed_thetas_eps(thetas_eps.n_elem); // changed
  arma::vec xs_eps = B_basis(xs, knots_t) * thetas_eps;
  arma::vec ws_eps(N);
  arma::vec us_eps(N);
  for (arma::uword NN = 0; NN < N; NN++) {
    us_eps(NN) = R::rnorm(0, 1);
    ws_eps(NN) = xs_eps(inds(NN)) + us_eps(NN);
  }


  arma::mat cond = arma::eye(K_t + 1, K_t + 1);
  for (arma::uword i = 1; i < K_t + 1; i++) {
    cond(i, i - 1) = -1;
  }
  arma::mat icond = arma::inv(cond);
  arma::vec cond_upr(K_t + 1, arma::fill::value(arma::datum::inf));
  arma::vec cond_lwr(K_t + 1, arma::fill::zeros);
  cond_lwr(0) = -arma::datum::inf;

    // ~~~~~~~~~~~~~~~~ //
  // For progress bar //
  // ~~~~~~~~~~~~~~~~ //
  arma::uword next_update = 5;

  // ~~~~~~~~~~~~~~~~~~~~~~~ //
  // Storage for MCMC output //
  // ~~~~~~~~~~~~~~~~~~~~~~~ //
  arma::mat ws_mcmc(nsamp, N);
  arma::vec current_xs;
  arma::vec proposed_xs;
  proposed_xs = current_xs = xs;
  arma::mat xs_mcmc(nsamp, n);
  arma::vec current_xs_eps;
  arma::vec proposed_xs_eps;
  proposed_xs_eps = current_xs_eps = xs_eps;
  arma::mat xs_eps_mcmc(nsamp, n);
  arma::vec xs_consumption_days(n);
  arma::vec current_xs_consumption_days;
  arma::vec proposed_xs_consumption_days;
  for(arma::uword nn = 0; nn<n; nn++){
  xs_consumption_days(nn) = xs(nn)/R::pnorm(xs_eps(nn), 0, 1, true, false);
  }

  arma::mat B_basis_store = B_basis(xs_grid, knots_t);
  arma::mat B_basis_var_grid_knots_t = B_basis(var_grid, knots_t);

  arma::vec current_thetas = thetas;
  arma::vec vars = B_basis(xs_consumption_days, knots_t) * exp(current_thetas);
  arma::vec current_vars = vars;

  proposed_xs_consumption_days = current_xs_consumption_days = xs_consumption_days;
  arma::mat xs_consumption_days_mcmc(nsamp, n);
  arma::vec current_us_eps;
  arma::vec proposed_us_eps;
  proposed_us_eps = current_us_eps = us_eps;
  arma::vec current_us;
  arma::vec proposed_us;
  proposed_us = current_us = us;
  arma::mat us_mcmc(nsamp, N);

  arma::vec xs_star(N);

  arma::vec density_es_est(es_grid.n_elem);
  arma::vec prob_consumption_est(x_grd_res);

  arma::mat current_likelihood_eps;
  arma::mat proposed_likelihood_eps;
  current_likelihood_eps = arma::mat(2, n, arma::fill::ones);
  proposed_likelihood_eps = current_likelihood_eps;
  arma::mat temp_current_us_likelihood_eps;
  arma::mat temp_proposed_us_likelihood_eps;
  temp_current_us_likelihood_eps = arma::mat(2, N, arma::fill::ones);
  temp_proposed_us_likelihood_eps = temp_current_us_likelihood_eps;

  arma::vec thetas_eps_est(thetas.n_elem);
  arma::vec thetas_est(thetas.n_elem);
  arma::mat thetas_xs_density_mcmc(nsamp, thetas_xs_density.n_elem,
                                   arma::fill::zeros);
  arma::mat thetas_mcmc(nsamp, thetas.n_elem, arma::fill::zeros);
  arma::mat thetas_eps_mcmc(nsamp, thetas_eps.n_elem,
                            arma::fill::zeros);

  double var_es = 1;
  arma::vec var_es_mcmc(nsamp);

  arma::vec current_thetas_xs_density;
  arma::vec proposed_thetas_xs_density;
  proposed_thetas_xs_density = current_thetas_xs_density = thetas_xs_density;

  arma::mat tempmat(n, x_grd_res);
  for (arma::uword rr = 0; rr < n; rr++) {
    for (arma::uword cc = 0; cc < x_grd_res; cc++) {
      tempmat(rr, cc) = abs(xs(rr) - xs_grid(cc));
    }
  }
  arma::uvec close_ind(n);
  for (arma::uword rr = 0; rr < n; rr++) {
    close_ind(rr) = tempmat.row(rr).index_min();
  }
  arma::vec d_ordinates_xs = B_basis_store.rows(close_ind) *
    B_basis_density_coeffs(thetas_xs_density, delta_t, K_t);
  arma::mat current_d_ordinates_xs;
  arma::mat proposed_d_ordinates_xs;
  proposed_d_ordinates_xs = current_d_ordinates_xs = d_ordinates_xs;

  arma::uvec z_us(N, arma::fill::ones);
  arma::vec pi_us(K_epsilon, arma::fill::value(1.0 / K_epsilon));
  arma::mat params_us(K_epsilon, 4);
  arma::rowvec params_us_ = {0.5, 0, 1,1};
  for (arma::uword kk = 0; kk < K_epsilon; kk++) {
    params_us.row(kk) = params_us_;
  }
  arma::umat z_us_mcmc(nsamp, N);
  arma::vec k_us_mcmc(nsamp);
  arma::mat pi_us_mcmc(nsamp, K_epsilon);
  arma::mat p_us_mcmc(nsamp, K_epsilon);
  arma::mat mu_us_mcmc(nsamp, K_epsilon);
  arma::mat sig1_us_mcmc(nsamp, K_epsilon);
  arma::mat sig2_us_mcmc(nsamp, K_epsilon);

  // ~~~~~~~~~~~~ //
  // Conduct MCMC //
  // ~~~~~~~~~~~~ //
  arma::uword niter = nsamp * nthin + nburn - nthin + 1;
  for (arma::uword it = 0; it < niter; it++) {
    // Updating thetas.xs.density
    arma::mat theta_xs_density_cov(K_t + 1, K_t + 1);
    theta_xs_density_cov.diag() = arma::vec(K_t + 1, arma::fill::value(sig_tune_thetas_1));
    theta_xs_density_cov += sig_tune_thetas_2 * prop_sig_thetas_xs_density;
    proposed_thetas_xs_density = rmvnorm(current_thetas_xs_density,
                                         theta_xs_density_cov).t(); // changed to prop_sig_thetas_xs_density from prop_sig_thetas

    for (arma::uword rr = 0; rr < n; rr++) {
      for (arma::uword cc = 0; cc < x_grd_res; cc++) {
        tempmat(rr, cc) = abs(xs(rr) - xs_grid(cc));
      }
    }
    for (arma::uword rr = 0; rr < n; rr++) {
      close_ind(rr) = tempmat.row(rr).index_min();
    }

    current_d_ordinates_xs = B_basis_store.rows(close_ind) *
      B_basis_density_coeffs(current_thetas_xs_density, delta_t, K_t);
    proposed_d_ordinates_xs = B_basis_store.rows(close_ind) *
      B_basis_density_coeffs(proposed_thetas_xs_density, delta_t, K_t);

    arma::vec current_log_prior = -current_thetas_xs_density.t() *
      P_t * current_thetas_xs_density / (2.0 * sigmasq_xi);
    arma::vec proposed_log_prior = -proposed_thetas_xs_density.t() *
      P_t * proposed_thetas_xs_density / (2.0 * sigmasq_xi);

    double current_log_likelihood = accu(log(current_d_ordinates_xs));
    double proposed_log_likelihood = accu(log(proposed_d_ordinates_xs));

    double log_mh_ratio = proposed_log_prior(0) + proposed_log_likelihood -
      current_log_likelihood - current_log_prior(0);
    if (log_mh_ratio == arma::datum::nan) {
      log_mh_ratio = -arma::datum::inf;
    }
    if (log(R::runif(0, 1)) < log_mh_ratio) {
      thetas_xs_density = current_thetas_xs_density =
        proposed_thetas_xs_density;
      d_ordinates_xs = current_d_ordinates_xs = proposed_d_ordinates_xs;
    }

    // Updating sigmasq_xi
    arma::vec sigmasq_xi_rate = b_xi + thetas_xs_density.t() *
      P_t * thetas_xs_density / 2.0;
    sigmasq_xi = 1.0 / R::rgamma(a_xi + (K_t + 1) / 2.0,
                                     1.0 / sigmasq_xi_rate(0));

    // Updating xs (and us)
    for (arma::uword nn = 0; nn < n; nn++) {
      arma::vec val = rtnorm(1, current_xs(nn), 0.1, xs_lwr, xs_upr);
      proposed_xs(nn) = val(0);
    }
    for (arma::uword rr = 0; rr < n; rr++) {
      for (arma::uword cc = 0; cc < x_grd_res; cc++) {
        tempmat(rr, cc) = abs(proposed_xs(rr) - xs_grid(cc));
      }
    }
    for (arma::uword rr = 0; rr < n; rr++) {
      close_ind(rr) = tempmat.row(rr).index_min();
    }
    arma::mat proposed_prior = B_basis_store.rows(close_ind) *
      B_basis_density_coeffs(thetas_xs_density, delta_t, K_t);
    arma::mat current_prior = d_ordinates_xs;

    arma::vec proposed_xs_eps = B_basis_store.rows(close_ind) * thetas_eps;
    for (arma::uword nn = 0; nn < n; nn++) {
      proposed_xs_consumption_days(nn) = proposed_xs(nn) /
        R::pnorm(proposed_xs_eps(nn), 0, 1, true, false);
    }
    arma::uvec excss_ = arma::find(proposed_xs_consumption_days > max(knots_t));
    proposed_xs_consumption_days(excss_) = arma::vec(excss_.n_elem,
                                 arma::fill::value(
                                   max(knots_t) - range_start_xs / 1000
                                 ));

    for (arma::uword rr = 0; rr < n; rr++) {
      for (arma::uword cc = 0; cc < x_grd_res; cc++) {
        tempmat(rr, cc) = abs(proposed_xs_consumption_days(rr) - xs_grid(cc));
      }
    }
    for (arma::uword rr = 0; rr < n; rr++) {
      close_ind(rr) = tempmat.row(rr).index_min();
    }
    arma::vec proposed_vars = B_basis_store.rows(close_ind) * exp(thetas);

    proposed_us_eps = ws_eps - proposed_xs_eps(inds);
    arma::vec temp_current_us_likelihood(N);
    arma::vec temp_proposed_us_likelihood(N);
    for (arma::uword NN = 0; NN < N; NN++) {
      temp_current_us_likelihood(NN) =
        R::dnorm(current_us_eps(NN), 0, 1, false);
      temp_proposed_us_likelihood(NN) =
        R::dnorm(proposed_us_eps(NN), 0, 1, false);
    }
    arma::vec current_likelihood_eps(n);
    arma::vec proposed_likelihood_eps(n);
    arma::uword m_sum = 0;
    for (arma::uword nn = 0; nn < n; nn++) {
      arma::uword m = mis(nn);
      current_likelihood_eps(nn) = prod(temp_current_us_likelihood(
        arma::span(m_sum, m_sum + m - 1)
      ));
      proposed_likelihood_eps(nn) = prod(temp_proposed_us_likelihood(
        arma::span(m_sum, m_sum + m - 1)
      ));
      m_sum += m;
    }

    proposed_us = ws - proposed_xs_consumption_days(inds);
    arma::uword k_us;
    if (type_u == 0) {
      k_us = max(z_us);
      temp_current_us_likelihood = fu_mixnorm(current_us, 0,
                                              sqrt(current_vars(inds)),
                                              pi_us(arma::span(0, k_us - 1)),
                                              params_us.rows(
                                                arma::span(0, k_us - 1)
                                              ));
      temp_proposed_us_likelihood = fu_mixnorm(proposed_us, 0,
                                              sqrt(proposed_vars(inds)),
                                              pi_us(arma::span(0, k_us - 1)),
                                              params_us.rows(
                                                arma::span(0, k_us - 1)
                                              ));
    } else if (type_u == 1) {
      for (arma::uword NN = 0; NN < N; NN++) {
        temp_current_us_likelihood(NN) = R::dnorm(current_us(NN), 0, sqrt(current_vars(inds(NN))), false);
        temp_proposed_us_likelihood(NN) = R::dnorm(proposed_us(NN), 0, sqrt(proposed_vars(inds(NN))), false);
      }
    } else if (type_u == 2) {
      for (arma::uword NN = 0; NN < N; NN++) {
        temp_current_us_likelihood(NN)  = 0.5 * R::dexp(abs(current_us(NN)), 1.0 / (sqrt(current_vars(inds(NN))) * 2.0), false);
        temp_proposed_us_likelihood(NN) = 0.5 * R::dexp(abs(proposed_us(NN)), 1.0 / (sqrt(proposed_vars(inds(NN))) * 2.0), false);
      }
    }

    arma::vec current_likelihood(n);
    arma::vec proposed_likelihood(n);
    m_sum = 0;
    for (arma::uword nn = 0; nn < n; nn++) {
      arma::uword m = mis(nn);
      current_likelihood(nn) = prod(temp_current_us_likelihood(
        arma::span(m_sum, m_sum + m - 1)
      ));
      proposed_likelihood(nn) = prod(temp_proposed_us_likelihood(
        arma::span(m_sum, m_sum + m - 1)
      ));
      m_sum += m;
    }
    if (true){
    current_likelihood %= current_likelihood_eps;
    proposed_likelihood %= proposed_likelihood_eps;
    }

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
    current_xs_eps(inds_to_replace) = proposed_xs_eps(inds_to_replace);
    xs_eps(inds_to_replace) = current_xs_eps(inds_to_replace);
    current_xs_consumption_days(inds_to_replace) =
      proposed_xs_consumption_days(inds_to_replace);
    xs_consumption_days(inds_to_replace) =
      current_xs_consumption_days(inds_to_replace);
    current_vars(inds_to_replace) = proposed_vars(inds_to_replace);
    vars(inds_to_replace) = current_vars(inds_to_replace);

    us_eps = current_us_eps = ws_eps - xs_eps(inds);
    us = current_us = ws - xs_consumption_days(inds);

    // Updating thetas
    if (type_v) {
      arma::mat theta_cov(K_t + 1, K_t + 1);
      theta_cov.diag() = arma::vec(K_t + 1, arma::fill::value(sig_tune_thetas_1));
      theta_cov += sig_tune_thetas_2 * prop_sig_thetas;
      arma::vec proposed_thetas = rmvnorm(current_thetas, theta_cov).t();
      proposed_thetas(0) =  2 * proposed_thetas(1) - proposed_thetas(2);
      for (arma::uword rr = 0; rr < n; rr++) {
        for (arma::uword cc = 0; cc < x_grd_res; cc++) {
          tempmat(rr, cc) = abs(xs_consumption_days(rr) - xs_grid(cc));
        }
      }
      for (arma::uword rr = 0; rr < n; rr++) {
        close_ind(rr) = tempmat.row(rr).index_min();
      }
      proposed_vars = B_basis_store.rows(close_ind) * exp(proposed_thetas);

      current_log_prior = -current_thetas.t() * P_t * current_thetas;
      current_log_prior /= (2.0 * sigmasq_vartheta);
      proposed_log_prior = - proposed_thetas.t() * P_t * proposed_thetas;
      proposed_log_prior /= (2.0 * sigmasq_vartheta);

      arma::vec mean_(N, arma::fill::zeros);
      arma::vec current_sds_rep(N);
      arma::vec proposed_sds_rep(N);
      arma::uword us_idx = 0;
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
          temp_current_us_likelihood(NN)  = 0.5 * R::dexp(abs(us(NN)), sqrt(2.0) / current_sds_rep(NN), false);
          temp_proposed_us_likelihood(NN) = 0.5 * R::dexp(abs(us(NN)), sqrt(2.0) / proposed_sds_rep(NN), false);
        }
      }

      current_log_likelihood = sum(log(temp_current_likelihood));
      proposed_log_likelihood = sum(log(temp_proposed_likelihood));

      log_mh_ratio = proposed_log_prior(0) + proposed_log_likelihood;
      log_mh_ratio -= current_log_likelihood + current_log_prior(0);
      if (log_mh_ratio == arma::datum::nan) {
        log_mh_ratio = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_ratio) {
        thetas = current_thetas = proposed_thetas;
        vars = current_vars = proposed_vars;
      }
    }

    // Update sigmasq_vartheta
    arma::vec sigmasq_vartheta_rate = b_vartheta + thetas.t() * P_t * thetas / 2.0;
    sigmasq_vartheta = 1.0 / R::rgamma(a_vartheta + (K_t + 1) / 2.0, 1.0 / sigmasq_vartheta_rate(0));

    if (((it >= 0) | (it >= (nburn / 2))) & (type_u == 0)) {
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
        arma::uvec samp_ = sample(K_epsilon, 1, false, prob_us);
        z_us(NN) = samp_(0);
      }

      // Update cluster probabilities
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


      k_us = max(z_us);
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
          params_us.row(kk) = r_proposal_params_restricted_mix_norm(1, 1, 3, 3,
                        3, 3, 3).t();
        }
      }

      var_es = var_e_fn(pi_us(arma::span(0, k_us - 1)),
                        params_us.rows(arma::span(0, k_us - 1)));
      params_us(arma::span(0, k_us - 1), 1) /= sqrt(var_es);
      params_us(arma::span(0, k_us - 1), arma::span(2, 3)) /= var_es;
    }

    // Update ws_eps
    xs_star = xs_eps(inds);
    for (arma::uword idx = 0; idx < inc.n_elem; idx++) {
      arma::vec val_ = rtnorm(1, xs_star(inc(idx)), 1, -arma::datum::inf, 0);
      ws_eps(inc(idx)) = val_(0);
    }
    for (arma::uword idx = 0; idx < ic.n_elem; idx++) {
      arma::vec val_ = rtnorm(1, xs_star(ic(idx)), 1, 0, arma::datum::inf);
      ws_eps(ic(idx)) = val_(0);
    }

    // Update thetas_eps
    if (true) {
      for (arma::uword rr = 0; rr < n; rr++) {
        for (arma::uword cc = 0; cc < x_grd_res; cc++) {
          tempmat(rr, cc) = abs(xs(rr) - xs_grid(cc));
        }
      }
      for (arma::uword rr = 0; rr < n; rr++) {
        close_ind(rr) = tempmat.row(rr).index_min();
      }
      arma::vec mu_thetas_eps(K_t + 1, arma::fill::zeros);
      arma::mat sigma_thetas_eps(K_t + 1, K_t + 1, arma::fill::zeros);
      for (arma::uword nn = 0; nn < n; nn++) {
        for (arma::uword jj = inds1(nn); jj <= inds2(nn); jj++) {
          mu_thetas_eps += ws_eps(jj) * B_basis_store.row(close_ind(nn)).t();
        }
        sigma_thetas_eps += mis(nn) * B_basis_store.row(close_ind(nn)).t() *
          B_basis_store.row(close_ind(nn));
      }
      sigma_thetas_eps = arma::inv(P_t / sigmasq_beta +
        arma::inv(sigma0_thetas_eps) + sigma_thetas_eps);
      mu_thetas_eps = sigma_thetas_eps *
        (arma::inv(sigma0_thetas_eps) * mu0_thetas_eps + mu_thetas_eps);
      arma::vec cond_mu = cond * mu_thetas_eps;
      proposed_thetas_eps = mu_thetas_eps +
        icond * rmvtnorm(cond_lwr - cond_mu,
                         cond_upr - cond_mu,
                         cond * sigma_thetas_eps * cond.t());


      proposed_xs_eps = B_basis_store.rows(close_ind) * proposed_thetas_eps;
      for (arma::uword nn = 0; nn < xs_consumption_days.n_elem; nn++) {
        proposed_xs_consumption_days(nn) = xs(nn) /
          R::pnorm(proposed_xs_eps(nn), 0, 1, true, false);
      }
      arma::vec mean_(N, arma::fill::zeros);
      arma::vec temp_current_likelihood(N);
      arma::vec temp_proposed_likelihood(N);
      proposed_us = ws - proposed_xs_consumption_days(inds);
      if (type_u == 0) {
        temp_current_likelihood = d_restricted_mixnorm(current_us, mean_,
                                                       sqrt(vars(inds)),
                                                       params_us.rows(
                                                         z_us - 1
                                                       ));
        temp_proposed_likelihood = d_restricted_mixnorm(proposed_us, mean_,
                                                        sqrt(vars(inds)),
                                                        params_us.rows(
                                                          z_us - 1
                                                        ));
      } else if (type_u == 1) {
        for (arma::uword NN = 0; NN < N; NN++) {
          temp_current_likelihood(NN) = R::dnorm(us(NN), 0, sqrt(vars(inds(NN))), false);
          temp_proposed_likelihood(NN) = R::dnorm(us(NN), 0, sqrt(vars(inds(NN))), false);
        }
      } else if (type_u == 2) {
        for (arma::uword NN = 0; NN < N; NN++) {
          temp_current_us_likelihood(NN)  = 0.5 * R::dexp(abs(us(NN)), sqrt(2.0) / sqrt(vars(inds(NN))), false);
          temp_proposed_us_likelihood(NN) = 0.5 * R::dexp(abs(us(NN)), sqrt(2.0) / sqrt(vars(inds(NN))), false);
        }
      }

      current_log_likelihood = sum(log(temp_current_likelihood));
      proposed_log_likelihood = sum(log(temp_proposed_likelihood));

      log_mh_ratio = proposed_log_likelihood - current_log_likelihood;
      if (log_mh_ratio == arma::datum::nan) {
        log_mh_ratio = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_ratio) {
        thetas_eps =  proposed_thetas_eps;
      }
      xs_eps = current_xs_eps = B_basis_store.rows(close_ind) * thetas_eps;
      for (arma::uword nn = 0; nn < xs_consumption_days.n_elem; nn++) {
        xs_consumption_days(nn) = current_xs_consumption_days(nn) = xs(nn) /
          R::pnorm(xs_eps(nn), 0, 1, true, false);
      }
      arma::uvec idx_ = arma::find(xs_consumption_days > max(knots_t));
      for (arma::uword nn = 0; nn < idx_.n_elem; nn++) {
        xs_consumption_days(idx_(nn)) = max(knots_t) - range_start_xs / 1000;
      }
      for (arma::uword rr = 0; rr < n; rr++) {
        for (arma::uword cc = 0; cc < x_grd_res; cc++) {
          tempmat(rr, cc) = abs(xs_consumption_days(rr) - xs_grid(cc));
        }
      }
      for (arma::uword rr = 0; rr < n; rr++) {
        close_ind(rr) = tempmat.row(rr).index_min();
      }
      vars = current_vars = B_basis_store.rows(close_ind) * exp(thetas);
      }

    // Update sigmasq_beta
    arma::vec sigmasq_beta_rate = b_beta + thetas_eps.t() * P_t * thetas_eps / 2.0;
    sigmasq_beta = 1.0 / R::rgamma(a_beta + (K_t + 1) / 2.0,
                              1.0 / sigmasq_beta_rate(0));

    // Update ws when not observed
    arma::vec mu_ws(N);
    arma::vec sigmasq_ws(N);
    if (type_u == 0) {
      arma::vec c2_us(N);
      for (arma::uword NN = 0; NN < N; NN++) {
        c2_us(NN) = R::rbinom(1, params_us(z_us(NN) - 1));
      }
      arma::mat mu_us(params_us.n_rows, 2);
      mu_us.col(0) = params_us.col(1) % (1.0 - params_us.col(0)) /
        (params_us.col(0) % params_us.col(0) +
          (1.0 - params_us.col(0)) % (1.0 - params_us.col(0)));
      mu_us.col(1) = - params_us.col(1) % params_us.col(0) /
        (params_us.col(0) % params_us.col(0) +
          (1.0 - params_us.col(0)) % (1.0 - params_us.col(0)));
      for (arma::uword NN = 0; NN < N; NN++) {
        mu_ws(NN) = xs_consumption_days(inds(NN));
        mu_ws(NN) += sqrt(vars(inds(NN))) * mu_us(z_us(NN) - 1, c2_us(NN));
        sigmasq_ws(NN) = vars(inds(NN)) * params_us(z_us(NN) - 1, c2_us(NN) + 2);
      }
    } else {
      for (arma::uword NN = 0; NN < N; NN++) {
        mu_ws(NN) = xs_consumption_days(inds(NN));
        sigmasq_ws(NN) = vars(inds(NN));
      }
    }
    for (arma::uword idx = 0; idx < inc.n_elem; idx++) {
      ws(inc(idx)) = R::rnorm(mu_ws(inc(idx)), sqrt(sigmasq_ws(inc(idx))));
    }

    us_eps = current_us_eps = ws_eps - xs_eps(inds);
    us = current_us = (ws - xs_consumption_days(inds));

    if ((it >= nburn) & ((it - nburn) % nthin == 0)) {
      arma::uword idx = (it - nburn) / nthin;
      ws_mcmc.row(idx) = ws.t();
      xs_mcmc.row(idx) = xs.t();
      xs_eps_mcmc.row(idx) = xs_eps.t();
      xs_consumption_days_mcmc.row(idx) = xs_consumption_days.t();
      us_mcmc.row(idx) = us.t();
      thetas_xs_density_mcmc.row(idx) = thetas_xs_density.t();
      thetas_mcmc.row(idx) = thetas.t();
      thetas_eps_mcmc.row(idx) = thetas_eps.t();

      if (type_u == 0) {
        var_es = var_e_fn(pi_us(arma::span(0, k_us - 1)),
                          params_us.rows(arma::span(0, k_us - 1)));
        z_us_mcmc.row(idx) = z_us.t();
        k_us_mcmc(idx) = max(z_us);
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

      thetas_est += log(var_es) + thetas;
      thetas_eps_est += thetas_eps;
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

  return Rcpp::List::create(Rcpp::_("w") = ws_mcmc,
                            Rcpp::_("x") = xs_mcmc,
                            Rcpp::_("theta_x") = thetas_xs_density_mcmc,
                            Rcpp::_("theta_v") = thetas_mcmc,
                            Rcpp::_("theta_e") = thetas_eps_mcmc,
                            Rcpp::_("var_e") = var_es_mcmc,
                            Rcpp::_("x_e") = xs_eps_mcmc,
                            Rcpp::_("x_nz") = xs_consumption_days_mcmc,
                            Rcpp::_("u") = us_mcmc,
                            Rcpp::_("k_u") = k_us_mcmc,
                            Rcpp::_("z_u") = z_us_mcmc,
                            Rcpp::_("pi_u") = pi_us_mcmc,
                            Rcpp::_("p_u") = p_us_mcmc,
                            Rcpp::_("mu_u") = mu_us_mcmc,
                            Rcpp::_("sig1_u") = sig1_us_mcmc,
                            Rcpp::_("sig2_u") = sig2_us_mcmc,
                            Rcpp::_("knots") = knots_t);
}
