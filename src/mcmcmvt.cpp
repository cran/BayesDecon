#include <math.h>
#include <RcppArmadillo.h>
#include "tapply.h"
#include "B_basis.h"
#include "tnorm.h"
#include "sample.h"
#include "rdirichlet.h"
#include "mvnorm.h"
#include "mixnorm.h"
#include "formcorr.h"
#include "var_e_fn.h"
#include "get_close.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List mcmcmvt_cpp(arma::uword nsamp, arma::uword nburn,
                       arma::uword nthin, arma::uword niter_u,
                       arma::mat ws, arma::mat xs, arma::mat us,
                       arma::uvec mis, arma::uword p, arma::uword q,
                       arma::mat knots_t, arma::vec delta_t,
                       arma::mat thetas, arma::mat theta_x,
                       arma::umat z_u, arma::mat pi_u,
                       arma::umat ic, arma::umat inc,
                       arma::mat x_e, arma::mat x_nz, arma::mat theta_e,
                       arma::umat z_x, arma::mat pi_x, arma::cube params_x,
                       arma::mat P_t, arma::cube prop_sig_thetas, arma::cube prop_sig_thetas_x_density,
                       arma::vec xs_lwr, arma::vec xs_upr,
                       arma::uword x_grd_res, arma::uword var_grd_res,
                       arma::uword e_grd_res,
                       double alpha_x,
                       double a_sigmasq_x, double b_sigmasq_x,
                       arma::uword K_x, double a_vartheta, double b_vartheta,
                       bool type_v, arma::uword K_t,
                       arma::uword type_u, arma::uword K_epsilon, double alpha_epsilon,
                       double a_epsilon, double b_epsilon, double sigmasq_epsilon,
                       arma::vec a_xi, arma::vec b_xi,
                       arma::vec sigmasq_xi,
                       arma::vec a_beta, arma::vec b_beta,
                       arma::vec sigmasq_beta,
                       arma::mat mu0_thetas_eps, double sigma00_thetas_eps,
                       double sig_tune_theta_1, double sig_tune_theta_2, bool show_progress) {
  // Set scalars
  arma::uword n = mis.n_elem;
  arma::uword N = sum(mis);
  arma::uword d = ws.n_rows;

  // Populate useful indices
  arma::uvec inds(N);
  arma::uvec inds1(n);
  arma::uvec inds2(n);
  inds1(0) = 1;
  inds2(0) = mis(0);
  arma::uword tmp_idx = 0;
  for (arma::uword nn = 0; nn < mis.n_elem; nn++) {
    if (nn > 0) {
      inds1(nn) = inds1(nn - 1) + mis(nn - 1);
      inds2(nn) = inds1(nn) + mis(nn) - 1;
    }
    for (arma::uword mm = 0; mm < mis(nn); mm++) {
      inds(tmp_idx) = nn;
      tmp_idx++;
    }
  }

  arma::uvec num_c = sum(ic, 0).t();
  arma::uvec num_nc = sum(inc, 0).t();

  // Condition Matrix for Monotone Coefficients
  arma::mat cond = arma::eye(K_t + 1, K_t + 1);
  for (arma::uword i = 1; i < K_t + 1; i++) {
    cond(i, i - 1) = -1;
  }
  arma::mat icond = arma::inv(cond);
  arma::vec cond_upr(K_t + 1, arma::fill::value(arma::datum::inf));
  arma::vec cond_lwr(K_t + 1, arma::fill::zeros);
  cond_lwr(0) = -arma::datum::inf;

  // Priors and Initial Values
  // // Initialize xs and us
  arma::mat s2is(d, n, arma::fill::zeros);
  arma::mat start_xs = xs;
  arma::mat curr_xs = xs;
  arma::mat prop_xs(d, n, arma::fill::zeros);
  arma::mat prop_x_nz(d, n, arma::fill::zeros);
  arma::mat curr_x_nz(d, n, arma::fill::zeros);
  arma::mat xs_trans(p, n, arma::fill::zeros);
  arma::mat x_grd(d, x_grd_res, arma::fill::zeros);
  arma::vec range_start_xs(d, arma::fill::zeros);
  for (arma::uword dd = 0; dd < d; dd++) {
    xs.row(dd).clamp(min(knots_t.row(dd)), max(knots_t.row(dd)));
    s2is.row(dd) = tapply_var(ws.row(dd).t(), mis).t();
    x_grd.row(dd) = arma::linspace<arma::rowvec>(xs_lwr(dd), xs_upr(dd),
               x_grd_res);
    range_start_xs(dd) = max(xs.row(dd)) - min(xs.row(dd));
  }

  double mu0_x = 0;
  double sigmasq0_x = 1;
  arma::umat n_kk_x(p, K_x, arma::fill::zeros);
  arma::mat prob_x(p, K_x, arma::fill::zeros);
  if (p > 0) {
    mu0_x = mean(mean(xs.rows(arma::span(q, q + p - 1)), 1));
    sigmasq0_x = mean(var(xs.rows(arma::span(q, q + p - 1)), 1, 1));
  }
  // // Prior and initialization of sigmasq_vartheta and thetas
  arma::vec sigmasq_vartheta(d, arma::fill::value(0.01));
  arma::vec var_es(d, arma::fill::ones);
  arma::mat curr_thetas = thetas;
  arma::mat prop_thetas(d, K_t + 1, arma::fill::zeros);
  arma::mat vars(d, n, arma::fill::zeros);
  arma::mat curr_vars(d, n, arma::fill::zeros);
  arma::mat prop_vars(d, n, arma::fill::zeros);
  arma::mat var_grd(d, var_grd_res, arma::fill::zeros);
  for (arma::uword dd = 0; dd < d; dd++) {
    var_grd.row(dd) = arma::linspace<arma::rowvec>(
      min(knots_t.row(dd)),
      max(knots_t.row(dd)),
      var_grd_res);
  }

  arma::umat indices_consumed = (ws != 0);
  arma::umat indices_not_consumed = (ws == 0);
  arma::uvec num_proxies_consumed = sum(indices_consumed, 0).t();
  arma::uvec num_proxies_not_consumed = sum(indices_not_consumed, 0).t();

  // Prior and initialization for mixture
  arma::umat n_kk_u(d, K_epsilon, arma::fill::zeros);
  arma::mat prob_u(d, K_epsilon, arma::fill::zeros);
  arma::vec e_grid = arma::linspace<arma::vec>(-4, 4, e_grd_res);
  arma::mat density_e_est(d, e_grd_res, arma::fill::zeros);

  //Initialize density parameters for episodic components
  arma::mat thetas_x_density = theta_x;
  arma::mat curr_thetas_x_density = thetas_x_density;
  arma::mat prop_thetas_x_density = thetas_x_density;

  // Initialize latent variables for episodic components
  arma::mat sigma0_thetas_eps = sigma00_thetas_eps *
    arma::eye(K_t + 1, K_t + 1);

  // Catch up after running univariate models
  arma::mat mu_x_mat = params_x(arma::span::all, arma::span::all, arma::span(0, 0));
  arma::mat sigmasq_x_mat = params_x(arma::span::all, arma::span::all,
                                     arma::span(1, 1));
  arma::vec mu_x(K_x);
  arma::vec sigmasq_x(K_x);
  for (arma::uword dd = 0; dd < d; dd++) {
    start_xs.row(dd) = xs.row(dd);
    curr_xs.row(dd) = start_xs.row(dd);
    arma::vec temp_vars = B_basis(x_nz.row(dd).t(), knots_t.row(dd).t()) *
      exp(thetas.row(dd)).t();
    temp_vars.replace(0, max(temp_vars));
    vars.row(dd) = temp_vars.t();
    curr_vars.row(dd) = vars.row(dd);
  }

  // Redefine shared component parameters mu_x and sigmasq_x
  if (p > 0) {
    double mu_x_min = arma::datum::inf;
    double mu_x_max = -arma::datum::inf;
    for (arma::uword pp = 0; pp < p; pp++) {
      arma::vec mu_x_ = mu_x_mat.row(pp).t();
      mu_x_min = std::min(mu_x_min, min(mu_x_(z_x.row(pp))));
      mu_x_max = std::max(mu_x_max, max(mu_x_(z_x.row(pp))));
    }
    arma::vec mu_x_new = arma::linspace<arma::vec>(mu_x_min, mu_x_max, K_x);
    for (arma::uword pp = 0; pp < p; pp++) {
      arma::mat tempmat(n, K_x);
      for (arma::uword rr = 0; rr < n; rr++) {
        for (arma::uword cc = 0; cc < K_x; cc++) {
          tempmat(rr, cc) = abs(mu_x_mat(pp, z_x(pp, rr)) - mu_x_new(cc));
        }
        z_x(pp, rr) = tempmat.row(rr).index_min();
      }
    }
    mu_x = mu_x_new;
    sigmasq_x = arma::vec(K_x,
                          arma::fill::value(2.0 *
                            var(arma::vectorise(xs)) / K_x
                          ));
    for (arma::uword it = 0; it < 100; it++) {
      for (arma::uword pp = 0; pp < p; pp++) {
        arma::vec xs_trans_tmp_num(n);
        arma::vec xs_trans_tmp_den(n);
        for (arma::uword nn = 0; nn < n; nn++) {
          arma::uword zz = z_x(pp, nn);
          double mu_ = mu_x(zz);
          double sigsq_ = sigmasq_x(zz);
          xs_trans_tmp_num(nn) = R::pnorm((xs(q + pp, nn) - mu_) / sqrt(sigsq_),
                           0, 1, true, false);
          xs_trans_tmp_num(nn) -= R::pnorm((xs_lwr(q + pp) - mu_) /
            sqrt(sigsq_), 0, 1, true, false);
          xs_trans_tmp_den(nn) = R::pnorm((xs_upr(q + pp) - mu_) / sqrt(sigsq_),
                           0, 1, true, false);
          xs_trans_tmp_den(nn) -= R::pnorm((xs_lwr(q + pp) - mu_) /
            sqrt(sigsq_), 0, 1, true, false);
          double xs_trans_tmp = xs_trans_tmp_num(nn) / xs_trans_tmp_den(nn);
          if (xs_trans_tmp >= 1) {
            xs_trans_tmp = 0.999;
          }
          if (xs_trans_tmp <= 0) {
            xs_trans_tmp = 0.001;
          }
          xs_trans_tmp = mu_ + sqrt(sigsq_) * R::qnorm(xs_trans_tmp, 0, 1, true,
                                    false);
          xs_trans(pp, nn) = xs_trans_tmp;
        }
      }
        for (arma::uword kk = 0; kk < K_x; kk++) {
          arma::vec xspool({});
          for (arma::uword pp = 0; pp < p; pp++) {
            arma::vec xs_trans_ = xs_trans.row(pp).t();
            arma::uvec tmp = arma::find(z_x.row(pp) == kk);
            xspool = arma::join_vert(xspool, xs_trans_(tmp));
          }
          // Update mu_x
          double post_sigmasq = 1.0 / (sum(n_kk_x.col(kk)) / sigmasq_x(kk) +
                                       1.0 / sigmasq0_x);
          double post_mu = (sum(xspool) / sigmasq_x(kk) + mu0_x / sigmasq0_x) *
            post_sigmasq;
          mu_x(kk) = R::rnorm(post_mu, sqrt(post_sigmasq));
          // Update sigmasq_x
          double post_a_sigsqx = a_sigmasq_x + xspool.n_elem / 2.0;
          double post_b_sigsqx = b_sigmasq_x +
            sum(pow(xspool - mu_x(kk), 2)) / 2.0;
          sigmasq_x(kk) = 1.0 / R::rgamma(post_a_sigsqx, 1.0 / post_b_sigsqx);
        }
      for (arma::uword pp = 0; pp < p; pp++) {
        for (arma::uword nn = 0; nn < n; nn++) {
          arma::vec prob_x = pi_x.row(pp).t();
          for (arma::uword kk = 0; kk < K_x; kk++) {
            prob_x(kk) *= dtnorm_sclr(xs(q + pp, nn),
                   mu_x(kk), sqrt(sigmasq_x(kk)),
                   xs_lwr(q + pp), xs_upr(q + pp), false);
          }
          if (arma::any(arma::vectorise(prob_x != prob_x || prob_x == arma::datum::inf || prob_x == -arma::datum::inf))){
            prob_x.replace(arma::datum::nan, 0);
            prob_x.replace(arma::datum::inf, 0);
          }
          arma::uvec val = sample(K_x, 1, true, prob_x);
          z_x(pp, nn) = val(0) - 1;
        }
        n_kk_x.row(pp) = arma::hist(z_x.row(pp).t(),
                   arma::linspace<arma::uvec>(0, K_x - 1, K_x)).t();
        pi_x.row(pp) = rdirichlet(
          arma::conv_to<arma::vec>::from(n_kk_x.row(pp).t()) + alpha_x / K_x
        ).t();
      }
    }
  }

  // Redefine shared component parameters params_u
  arma::mat params_u(K_epsilon, 4);
  arma::rowvec params_u_ = {0.5, 0, 1, 1};
  for (arma::uword kk = 0; kk < K_epsilon; kk++) {
    params_u.row(kk) = params_u_;
  }
  if (type_u == 0) {
    for (arma::uword it = 0; it < 100; it++) {
        for (arma::uword rr = 0; rr < niter_u; rr++) {
          for (arma::uword kk = 0; kk < K_epsilon; kk++) {
            arma::vec uspool({});
            arma::vec varspool({});
            for (arma::uword ddd = 0; ddd < d; ddd++) {
              arma::vec u_ = us.row(ddd).t();
              arma::vec var_ = vars.row(ddd).t();
              arma::uvec tmp = arma::find(z_u.row(ddd) == kk);
              uspool = arma::join_vert(uspool, u_(tmp));
              varspool = arma::join_vert(varspool, var_(inds(tmp)));
            }
            arma::vec prop_params_u =
              r_tnorm_proposal_params_restricted_mix_norm(params_u.row(kk).t());
            double prop_log_prior = R::dnorm(prop_params_u(1),
                                            0, sqrt(sigmasq_epsilon), true) -
                                              (a_epsilon + 1) *
                                              (log(prop_params_u(2)) +
                                                log(prop_params_u(3))) -
                                              b_epsilon *
                                              (1.0 / prop_params_u(2) +
                                                1.0 / prop_params_u(3));
            double curr_log_prior = R::dnorm(params_u(kk, 1),
                                            0, sqrt(sigmasq_epsilon), true) -
                                              (a_epsilon + 1) *
                                              (log(params_u(kk, 2)) +
                                                log(params_u(kk, 3))) -
                                              b_epsilon *
                                              (1.0 / params_u(kk, 2) +
                                                1.0 / params_u(kk, 3));
            arma::vec mu_(uspool.n_elem, arma::fill::zeros);
            arma::vec tmp_prop_loglik = log(d_restricted_mixnorm(
              uspool, mu_, sqrt(varspool), prop_params_u.t()
            ));
            arma::vec tmp_curr_loglik = log(d_restricted_mixnorm(
              uspool, mu_, sqrt(varspool), params_u.row(kk)
            ));
            double prop_loglik = sum(tmp_prop_loglik);
            double curr_loglik = sum(tmp_curr_loglik);
            double log_acc_prob = prop_log_prior + prop_loglik -
              curr_log_prior - curr_loglik;
            if (log(R::runif(0, 1)) < log_acc_prob) {
              params_u.row(kk) = prop_params_u.t();
            }
          }
        }
        for (arma::uword dd = 0; dd < d; dd++) {
        for (arma::uword NN = 0; NN < N; NN++) {
          arma::vec prob_u = pi_u.row(dd).t();
          prob_u %= d_restricted_mixnorm_sclr(us(dd, NN), 0,
                                              sqrt(vars(dd, inds(NN))), params_u);
          if (arma::any(arma::vectorise(prob_u != prob_u || prob_u == arma::datum::inf || prob_u == -arma::datum::inf))){
            prob_u.replace(arma::datum::nan, 0);
            prob_u.replace(arma::datum::inf, 0);
          }
          arma::uvec val = sample(K_epsilon, 1, true, prob_u);
          z_u(dd, NN) = val(0) - 1;
        }
        n_kk_u.row(dd) = arma::hist(z_u.row(dd).t(),
                  arma::linspace<arma::uvec>(0, K_epsilon - 1, K_epsilon)).t();
        pi_u.row(dd) = rdirichlet(
          arma::conv_to<arma::vec>::from(n_kk_u.row(dd).t()) + alpha_epsilon / K_epsilon
        ).t();
      }
    }
  }

  // Precompute and store B-spline bases
  arma::cube B_basis_var_grid_knots_t(d, var_grd_res, K_t + 1);
  arma::cube B_basis_store(d, x_grd_res, K_t + 1);
  for (arma::uword dd = 0; dd < d; dd++) {
    arma::mat tmp = B_basis(var_grd.row(dd).t(),
                            knots_t.row(dd).t());
    B_basis_var_grid_knots_t.row(dd) = B_basis(var_grd.row(dd).t(),
                                 knots_t.row(dd).t());
    B_basis_store.row(dd) = B_basis(x_grd.row(dd).t(), knots_t.row(dd).t());
  }

  // Initialize x_e, w_e, and u_e
  arma::mat w_e(q, N, arma::fill::zeros);
  arma::mat u_e(q, N, arma::fill::zeros);
  arma::mat curr_u_e;
  arma::mat prop_u_e;
  if (q > 0) {
    for (arma::uword qq = 0; qq < q; qq++) {
      arma::vec x_e_ = x_e.row(qq).t();
      arma::vec x_star = x_e_(inds);
      arma::uvec inc_tmp = arma::find(inc.row(qq).t() == 1);
      arma::uvec ic_tmp = arma::find(ic.row(qq).t() == 1);
      for (arma::uword NN = 0; NN < N; NN++) {
        arma::vec val;
        if (inc(qq, NN) == 1) {
          val = rtnorm(1, x_star(NN), 1, -arma::datum::inf, 0);
        } else {
          val = rtnorm(1, x_star(NN), 1, 0, arma::datum::inf);
        }
        w_e(qq, NN) = val(0);
      }
    }
    u_e = w_e - x_e.submat(arma::linspace<arma::uvec>(0, q - 1, q), inds);
  }
  curr_u_e = u_e;
  prop_u_e = u_e;

  // Initialize copula and correlation matrix
  arma::mat yxs(d, n, arma::fill::zeros);
  arma::mat curr_yxs;
  arma::mat prop_yxs;
  if (q > 0) {
    for (arma::uword qq = 0; qq < q; qq++) {
      arma::vec tmp_xsrs = dist_norm_basis(xs.row(qq).t(), knots_t.row(qq).t(),
                                            thetas_x_density.row(qq).t());
      tmp_xsrs.clamp(0.001, 0.999);
      for (arma::uword nn = 0; nn < n; nn++) {
        yxs(qq, nn) = R::qnorm(tmp_xsrs(nn), 0, 1, true, false);
      }
    }
  }
  if (p > 0) {
    for (arma::uword pp = 0; pp < p; pp++) {
      arma::vec tmp_xsrs = dist_mixtnorm(xs.row(q + pp).t(), pi_x.row(pp).t(),
                                          mu_x, sqrt(sigmasq_x),
                                          xs_lwr(q + pp), xs_upr(q + pp));
      tmp_xsrs.clamp(0.001, 0.999);
      for (arma::uword nn = 0; nn < n; nn++) {
        yxs(q + pp, nn) = R::qnorm(tmp_xsrs(nn), 0, 1, true, false);
      }
    }
  }
  curr_yxs = yxs;
  prop_yxs = yxs;

  arma::mat yus(d, N, arma::fill::zeros);
  arma::mat curr_yus;
  arma::mat prop_yus;
  for (arma::uword dd = 0; dd < d; dd++) {
    arma::vec tmp_usrs(N);
    if (type_u == 0) {
      arma::uvec dd_(1, arma::fill::value(dd));
      tmp_usrs = dist_mixnorm(us.row(dd).t(),
                              0, sqrt(vars.submat(dd_, inds).t()),
                              pi_u.row(dd).t(), params_u);
    } else if (type_u == 1) {
      for (arma::uword NN = 0; NN < N; NN++) {
        tmp_usrs(NN) = R::dnorm(us(dd, NN), 0, sqrt(vars(dd, inds(NN))), false);
      }
    } else if (type_u == 2) {
      for (arma::uword NN = 0; NN < N; NN++) {
        tmp_usrs(NN) = R::dexp(abs(us(dd, NN)), 1.0 / sqrt(vars(dd, inds(NN))) / 2.0, false);
      }
    }
    tmp_usrs.clamp(0.001, 0.999);
    for (arma::uword NN = 0; NN < N; NN++) {
      yus(dd, NN) = R::qnorm(tmp_usrs(NN), 0, 1, true, false);
    }
  }
  curr_yus = yus;
  prop_yus = yus;

  arma::vec bx(d - 1, arma::fill::zeros);
  arma::vec curr_bx = bx;
  arma::vec thetax((d - 1) * (d - 2) / 2, arma::fill::zeros);
  arma::vec curr_thetax = thetax;
  arma::vec bxset = arma::linspace<arma::vec>(-0.99, 0.99, 41);
  arma::vec thetaxset = arma::linspace<arma::vec>(-3.14, 3.14, 41);

  arma::mat corr_x = formcorr(curr_bx, curr_thetax);
  arma::mat icorr_x = arma::inv(corr_x);
  arma::mat curr_corr_x = corr_x;
  arma::mat curr_icorr_x = icorr_x;

  arma::vec be(d - 1, arma::fill::zeros);
  arma::vec curr_be = be;
  arma::vec thetae((d - 1) * (d - 2) / 2, arma::fill::zeros);
  arma::vec curr_thetae = thetae;
  arma::vec beset = arma::linspace<arma::vec>(-0.99, 0.99, 41);
  arma::vec thetaeset = arma::linspace<arma::vec>(-3.14, 3.14, 41);

  arma::mat corr_e = formcorr(curr_be, curr_thetae);
  arma::mat icorr_e = arma::inv(corr_e);
  arma::mat curr_corr_e = corr_e;
  arma::mat curr_icorr_e = icorr_e;

  arma::mat d_ord_x(q, n);
  arma::mat curr_d_ord_x = d_ord_x;
  arma::mat prop_d_ord_x = d_ord_x;

  arma::mat curr_prior(d, n, arma::fill::ones);
  arma::mat prop_prior(d, n, arma::fill::ones);

  arma::mat curr_x_e = x_e;
  arma::mat prop_x_e = x_e;

  arma::mat prop_pc(d, n, arma::fill::zeros);
  arma::mat prop_cp(d, n, arma::fill::zeros);

  arma::mat prop_us = us;
  arma::mat curr_us = us;

  arma::mat curr_marglik(d, n, arma::fill::ones);
  arma::mat prop_marglik(d, n, arma::fill::ones);
  arma::mat curr_marglik_e(d, n, arma::fill::ones);
  arma::mat prop_marglik_e(d, n, arma::fill::ones);

  arma::vec curr_cop_prior(n, arma::fill::ones);
  arma::vec prop_cop_prior(n, arma::fill::ones);

  arma::vec curr_cop_like(n, arma::fill::ones);
  arma::vec prop_cop_like(n, arma::fill::ones);

  // ~~~~~~~~~~~~~~~~ //
  // For progress bar //
  // ~~~~~~~~~~~~~~~~ //
  arma::uword next_update = 5;

  //~~~~~~~~~~~~~~~~~~~~~~~~~//
  // Storage for MCMC Output //
  //~~~~~~~~~~~~~~~~~~~~~~~~~//
  arma::cube ws_mcmc(nsamp, d, N);
  arma::cube xs_mcmc(nsamp, d, n);
  arma::ucube z_x_mcmc(nsamp, p, n);
  arma::cube pi_x_mcmc(nsamp, p, K_x);
  arma::cube mu_x_mcmc(nsamp, p, K_x);
  arma::cube sig_x_mcmc(nsamp, p, K_x);
  arma::cube x_e_mcmc(nsamp, q, n);
  arma::cube x_nz_mcmc(nsamp, d, n);
  arma::cube thetas_x_dens_mcmc(nsamp, q, K_t + 1);
  arma::cube thetas_mcmc(nsamp, d, K_t + 1);
  arma::cube theta_e_mcmc(nsamp, d, K_t + 1);
  arma::mat var_es_mcmc(nsamp, d);
  arma::cube us_mcmc(nsamp, d, N);
  arma::umat k_u_mcmc(nsamp, d);
  arma::ucube z_u_mcmc(nsamp, d, N);
  arma::cube pi_u_mcmc(nsamp, d, K_epsilon);
  arma::mat p_u_mcmc(nsamp, K_epsilon);
  arma::mat mu_u_mcmc(nsamp, K_epsilon);
  arma::mat sig1_u_mcmc(nsamp, K_epsilon);
  arma::mat sig2_u_mcmc(nsamp, K_epsilon);

  arma::cube corr_x_mcmc(nsamp, d, d);
  arma::cube corr_e_mcmc(nsamp, d, d);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

  // Start MCMC
  arma::uword niter = nsamp * nthin + nburn - nthin + 1;
  for (arma::uword it = 0; it < niter; it++) {
    if (q > 0) {
      for (arma::uword qq = 0; qq < q; qq++) {
        // Update thetas_xs_density
        arma::mat theta_cov(K_t + 1, K_t + 1);
        theta_cov.diag() = arma::vec(K_t + 1,
                       arma::fill::value(sig_tune_theta_1));
        theta_cov += sig_tune_theta_2 * prop_sig_thetas.row(qq);
        prop_thetas_x_density.row(qq) = rmvnorm(
          curr_thetas_x_density.row(qq).t(),
          theta_cov
          );

        arma::uvec close_ind = get_close(xs.row(qq), x_grd.row(qq));

        arma::mat B_basis_store_qq = B_basis_store.row(qq);
        curr_d_ord_x.row(qq) = (B_basis_store_qq.rows(close_ind) *
          B_basis_density_coeffs(curr_thetas_x_density.row(qq).t(),
                                 delta_t(qq), K_t)).t();
        prop_d_ord_x.row(qq) = (B_basis_store_qq.rows(close_ind) *
          B_basis_density_coeffs(prop_thetas_x_density.row(qq).t(),
                                 delta_t(qq), K_t)).t();

        arma::vec curr_log_prior = -curr_thetas_x_density.row(qq) * P_t *
          curr_thetas_x_density.row(qq).t() / (2.0 * sigmasq_xi(qq));
        arma::vec prop_log_prior = -prop_thetas_x_density.row(qq) * P_t *
          prop_thetas_x_density.row(qq).t() / (2.0 * sigmasq_xi(qq));

        double curr_log_lik = sum(log(curr_d_ord_x.row(qq)));
        double prop_log_lik = sum(log(prop_d_ord_x.row(qq)));

        double log_mh_rat = prop_log_prior(0) + prop_log_lik -
          curr_log_lik - curr_log_prior(0);
        if (log_mh_rat == arma::datum::nan) {
          log_mh_rat = -arma::datum::inf;
        }
        if (log(R::runif(0, 1)) < log_mh_rat) {
          curr_thetas_x_density.row(qq) = prop_thetas_x_density.row(qq);
          thetas_x_density.row(qq) = curr_thetas_x_density.row(qq);
          curr_d_ord_x.row(qq) = prop_d_ord_x.row(qq);
          d_ord_x.row(qq) = curr_d_ord_x.row(qq);

          // Update yxs
          arma::vec tmp_xsrs = dist_norm_basis(xs.row(qq).t(),
                                               knots_t.row(qq).t(),
                                               thetas_x_density.row(qq).t());
          tmp_xsrs.clamp(0.001, 0.999);
          for (arma::uword nn = 0; nn < n; nn++) {
            curr_yxs(qq, nn) = R::qnorm(tmp_xsrs(nn), 0, 1, true, false);
            yxs(qq, nn) = curr_yxs(qq, nn);
          }
        }
        // Update s2t_x_density
        arma::vec sigmasq_xi_rate = b_xi(qq) + thetas_x_density.row(qq) *
          P_t * thetas_x_density.row(qq).t() / 2.0;
        sigmasq_xi(qq) = 1.0 / R::rgamma(a_xi(qq) + (K_t + 1) / 2.0,
                   1.0 / sigmasq_xi_rate(0));
      }
    }
    if (p > 0) {
      // Update mu_x and sigmasq_x
      for (arma::uword kk = 0; kk < K_x; kk++) {
        arma::vec xspool({});
        for (arma::uword pp = 0; pp < p; pp++) {
          arma::uvec idx = arma::find(z_x.row(pp) == kk);
          arma::vec xs_trans_ = xs_trans.row(pp).t();
          xspool = arma::join_vert(xspool, xs_trans_(idx));
        }

        double sigmasq_post = 1.0 / (sum(n_kk_x.col(kk)) / sigmasq_x(kk) +
                                    1.0 / sigmasq0_x);
        double mu_post = (sum(xspool) / sigmasq_x(kk) + mu0_x / sigmasq0_x) *
          sigmasq_post;
        mu_x(kk) = R::rnorm(mu_post, sqrt(sigmasq_post));

        double post_a_sigsqx = a_sigmasq_x + xspool.n_elem / 2.0;
        double post_b_sigsqx = b_sigmasq_x +
          sum(pow(xspool - mu_x(kk), 2)) / 2.0;
        sigmasq_x(kk) = 1.0 / R::rgamma(post_a_sigsqx, 1.0 / post_b_sigsqx);
      }
      // Update z_x and pi_x
      for (arma::uword pp = 0; pp < p; pp++) {
        for (arma::uword nn = 0; nn < n; nn++) {
          arma::vec prob_x = pi_x.row(pp).t();
          for (arma::uword kk = 0; kk < K_x; kk++) {
            prob_x(kk) *= dtnorm_sclr(xs(q + pp, nn),
                   mu_x(kk), sqrt(sigmasq_x(kk)),
                   xs_lwr(q + pp), xs_upr(q + pp), false);
          }
          if (arma::any(arma::vectorise(prob_x != prob_x || prob_x == arma::datum::inf || prob_x == -arma::datum::inf))){
            prob_x.replace(arma::datum::nan, 0);
            prob_x.replace(arma::datum::inf, 0);
          }
          arma::uvec val = sample(K_x, 1, false, prob_x);
          z_x(pp, nn) = val(0) - 1;
        }
        n_kk_x.row(pp) = arma::hist(z_x.row(pp).t(),
                   arma::linspace<arma::uvec>(0, K_x - 1, K_x)).t();
        pi_x.row(pp) = rdirichlet(
          arma::conv_to<arma::vec>::from(n_kk_x.row(pp).t()) + alpha_x / K_x
        ).t();

        // Update yxs
        arma::vec tmp_xsrs = dist_mixtnorm(xs.row(q + pp).t(),
                                           pi_x.row(pp).t(),
                                           mu_x,
                                           sqrt(sigmasq_x),
                                           xs_lwr(q + pp),
                                           xs_upr(q + pp));
        tmp_xsrs.clamp(0.001, 0.999);
        for (arma::uword nn = 0; nn < n; nn++) {
          curr_yxs(q + pp, nn) = R::qnorm(tmp_xsrs(nn), 0, 1, true, false);
          yxs(q + pp, nn) = curr_yxs(q + pp, nn);
        }
      }
    }



    // Update the Gaussian copula correlation parameters
    if (true) {//it >= 0
    for (arma::uword tt = 0; tt < (d - 1); tt++) {
      arma::vec prop_bx = curr_bx;
      prop_bx(tt) = sample(arma::vec({curr_bx(tt),
                           curr_bx(tt) + -2 * 0.99 / 40,
                           curr_bx(tt) + 2 * 0.99 / 40}));
      while(abs(prop_bx(tt)) > 0.99) {
        prop_bx(tt) = sample(arma::vec({curr_bx(tt),
                             curr_bx(tt) + -2 * 0.99 / 40,
                             curr_bx(tt) + 2 * 0.99 / 40}));
      }
      arma::mat prop_corr_x = formcorr(prop_bx, thetax);
      arma::mat prop_icorr_x = arma::inv(prop_corr_x);
      double prop_loglik = -(0.5 * n) * log(1 - prop_bx(tt) * prop_bx(tt));
      double curr_loglik = -(0.5 * n) * log(1 - curr_bx(tt) * curr_bx(tt));
      for (arma::uword nn = 0; nn < n; nn++) {
        arma::vec prop_val = yxs.col(nn).t() * prop_icorr_x * yxs.col(nn);
        arma::vec curr_val = yxs.col(nn).t() * curr_icorr_x * yxs.col(nn);
        prop_loglik -= 0.5 * prop_val(0);
        curr_loglik -= 0.5 * curr_val(0);
      }
      double log_mh_rat = prop_loglik - curr_loglik;
      if (log_mh_rat == arma::datum::nan) {
        log_mh_rat = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_rat) {
        bx(tt) = curr_bx(tt) = prop_bx(tt);
        corr_x = curr_corr_x = prop_corr_x;
        icorr_x = curr_icorr_x = prop_icorr_x;
      }
    }

    for (arma::uword tt = 0; tt < (d - 1) * (d - 2) / 2; tt++) {
      arma::vec prop_thetax = curr_thetax;
      prop_thetax(tt) = sample(arma::vec({curr_thetax(tt),
                               curr_thetax(tt) + -2 * 3.14 / 40,
                               curr_thetax(tt) + 2 * 3.14 / 40}));
      while (abs(prop_thetax(tt)) > 3.14) {
          prop_thetax(tt) = sample(arma::vec({curr_thetax(tt),
                                 curr_thetax(tt) + -2 * 3.14 / 40,
                                 curr_thetax(tt) + 2 * 3.14 / 40}));
      }
      arma::mat prop_corr_x = formcorr(bx, prop_thetax);
      arma::mat prop_icorr_x = arma::inv(prop_corr_x);

      double prop_loglik = 0;
      double curr_loglik = 0;
      for (arma::uword nn = 0; nn < n; nn++) {
        arma::vec prop_val = yxs.col(nn).t() * prop_icorr_x * yxs.col(nn);
        arma::vec curr_val = yxs.col(nn).t() * curr_icorr_x * yxs.col(nn);
        prop_loglik -= 0.5 * prop_val(0);
        curr_loglik -= 0.5 * curr_val(0);
      }
      double log_mh_rat = prop_loglik - curr_loglik;
      if (log_mh_rat == arma::datum::nan) {
        log_mh_rat = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_rat) {
        thetax(tt) = curr_thetax(tt) = prop_thetax(tt);
        corr_x = curr_corr_x = prop_corr_x;
        icorr_x = curr_icorr_x = prop_icorr_x;
      }
    }

    for (arma::uword tt = 0; tt < (d - 1); tt++) {
      arma::vec prop_be = curr_be;
      prop_be(tt) = sample(arma::vec({curr_be(tt),
                           curr_be(tt) + -2 * 0.99 / 40,
                           curr_be(tt) + 2 * 0.99 / 40}));
      while(abs(prop_be(tt)) > 0.99) {
        prop_be(tt) = sample(arma::vec({curr_be(tt),
                             curr_be(tt) + -2 * 0.99 / 40,
                             curr_be(tt) + 2 * 0.99 / 40}));
      }
      arma::mat prop_corr_e = formcorr(prop_be, thetae);
      arma::mat prop_icorr_e = arma::inv(prop_corr_e);

      double prop_loglik = -(0.5 * N) * log(1 - prop_be(tt) * prop_be(tt));
      double curr_loglik = -(0.5 * N) * log(1 - curr_be(tt) * curr_be(tt));
      for (arma::uword NN = 0; NN < N; NN++) {
        arma::vec prop_val = yus.col(NN).t() * prop_icorr_e * yus.col(NN);
        arma::vec curr_val = yus.col(NN).t() * curr_icorr_e * yus.col(NN);
        prop_loglik -= 0.5 * prop_val(0);
        curr_loglik -= 0.5 * curr_val(0);
      }
      double log_mh_rat = prop_loglik - curr_loglik;
      if (log_mh_rat == arma::datum::nan) {
        log_mh_rat = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_rat) {
        be(tt) = curr_be(tt) = prop_be(tt);
        corr_e = curr_corr_e = prop_corr_e;
        icorr_e = curr_icorr_e = prop_icorr_e;
      }
    }

    for (arma::uword tt = 0; tt < (d - 1) * (d - 2) / 2; tt++) {
      arma::vec prop_thetae = curr_thetae;
      prop_thetae(tt) = sample(arma::vec({curr_thetae(tt),
                               curr_thetae(tt) + -2 * 3.14 / 40,
                               curr_thetae(tt) + 2 * 3.14 / 40}));
      while (abs(prop_thetae(tt)) > 3.14) {
        prop_thetae(tt) = sample(arma::vec({curr_thetae(tt),
                                 curr_thetae(tt) + -2 * 3.14 / 40,
                                 curr_thetae(tt) + 2 * 3.14 / 40}));
      }
      arma::mat prop_corr_e = formcorr(be, prop_thetae);
      arma::mat prop_icorr_e = arma::inv(prop_corr_e);

      double prop_loglik = 0;
      double curr_loglik = 0;
      for (arma::uword NN = 0; NN < N; NN++) {
        arma::vec prop_val = yus.col(NN).t() * prop_icorr_e * yus.col(NN);
        arma::vec curr_val = yus.col(NN).t() * curr_icorr_e * yus.col(NN);
        prop_loglik -= 0.5 * prop_val(0);
        curr_loglik -= 0.5 * curr_val(0);
      }
      double log_mh_rat = prop_loglik - curr_loglik;
      if (log_mh_rat == arma::datum::nan) {
        log_mh_rat = -arma::datum::inf;
      }
      if (log(R::runif(0, 1)) < log_mh_rat) {
        thetae(tt) = curr_thetae(tt) = prop_thetae(tt);
        corr_e = curr_corr_e = prop_corr_e;
        icorr_e = curr_icorr_e = prop_icorr_e;
      }
    }
    }

    if (it >= nburn/2) {
      // Update xs and us
      if (q > 0) {
        for (arma::uword qq = 0; qq < q; qq++) {
          prop_xs.row(qq) = rtnorm(curr_xs.row(qq).t(), 0.1,
                      xs_lwr(qq), xs_upr(qq)).t(); // changed 0.1 to. sd_
          arma::mat tempmat(n, x_grd_res);
          arma::uvec close_ind(n);
          for (arma::uword rr = 0; rr < n; rr++) {
            for (arma::uword cc = 0; cc < x_grd_res; cc++) {
              tempmat(rr, cc) = abs(prop_xs(qq, rr) - x_grd(qq, cc));
            }
          }
          for (arma::uword rr = 0; rr < n; rr++) {
            close_ind(rr) = tempmat.row(rr).index_min();
          }

          arma::mat B_basis_store_qq = B_basis_store.row(qq);
          prop_prior.row(qq) = (B_basis_store_qq.rows(close_ind) *
            B_basis_density_coeffs(thetas_x_density.row(qq).t(),
                                   delta_t(qq),
                                   K_t)).t();
          curr_prior.row(qq) = d_ord_x.row(qq);

          prop_x_e.row(qq) = (B_basis_store_qq.rows(close_ind) *
            theta_e.row(qq).t()).t();
          prop_x_nz.row(qq) = prop_xs.row(qq);
          for (arma::uword nn = 0; nn < n; nn++) {
            prop_x_nz(qq, nn) /= R::pnorm(prop_x_e(qq, nn), 0, 1, true, false);
          }
           prop_x_nz.row(qq).clamp(0, max(knots_t.row(qq)) - 0.001);

          for (arma::uword rr = 0; rr < n; rr++) {
            for (arma::uword cc = 0; cc < x_grd_res; cc++) {
              tempmat(rr, cc) = abs(prop_x_nz(qq, rr) - x_grd(qq, cc));
            }
          }
          for (arma::uword rr = 0; rr < n; rr++) {
            close_ind(rr) = tempmat.row(rr).index_min();
          }
          prop_vars.row(qq) = (B_basis_store_qq.rows(close_ind) *
            exp(thetas.row(qq).t())).t();
          prop_vars.row(qq).replace(0, max(prop_vars.row(qq))); //why??

           double sd_ = (max(start_xs.row(qq)) - min(start_xs.row(qq))) / 6; // changed
          prop_pc.row(qq) = dtnorm(curr_xs.row(qq).t(), prop_xs.row(qq).t(),
                      sd_, xs_lwr(qq), xs_upr(qq), false).t();
          prop_cp.row(qq) = dtnorm(prop_xs.row(qq).t(), curr_xs.row(qq).t(),
                      sd_, xs_lwr(qq), xs_upr(qq), false).t();

          arma::vec tmp_xsrs = dist_norm_basis(xs.row(qq).t(),
                                               knots_t.row(qq).t(),
                                               thetas_x_density.row(qq).t());
          tmp_xsrs.clamp(0.001, 0.999);
          for (arma::uword nn = 0; nn < n; nn++) {
            yxs(qq, nn) = R::qnorm(tmp_xsrs(nn), 0, 1, true, false);
          }

          arma::vec prop_yxs_tmp = dist_norm_basis(
            prop_xs.row(qq).t(),
            knots_t.row(qq).t(),
            thetas_x_density.row(qq).t()
            );
          prop_yxs_tmp.clamp(0.001, 0.999);
          for (arma::uword nn = 0; nn < n; nn++) {
            prop_yxs(qq, nn) = R::qnorm(prop_yxs_tmp(nn), 0, 1, true, false);
          }

          arma::vec prop_xnz_qq = prop_x_nz.row(qq).t();
          prop_us.row(qq) = ws.row(qq) - prop_xnz_qq(inds).t();

          arma::vec prop_x_e_qq = prop_x_e.row(qq).t();
          prop_u_e.row(qq) = w_e.row(qq) - prop_x_e_qq(inds).t();
          arma::uword NN = 0;
          for (arma::uword nn = 0; nn < n; nn++) {
            double curr_tmp = 1;
            double prop_tmp = 1;
            for (arma::uword mm = 0; mm < mis(nn); mm++) {
              curr_tmp *= R::dnorm(curr_u_e(qq, NN), 0, 1, false);
              prop_tmp *= R::dnorm(prop_u_e(qq, NN), 0, 1, false);
              NN++;
            }
            curr_marglik_e(qq, nn) = curr_tmp;
            prop_marglik_e(qq, nn) = prop_tmp;
          }
        }
      }
      if (p > 0) {
        for (arma::uword pp = 0; pp < p; pp++) {
          double sd_ = (max(start_xs.row(q + pp)) -
                        min(start_xs.row(q + pp))) / 6;
          prop_xs.row(q + pp) = rtnorm(curr_xs.row(q + pp).t(), sd_,
                      xs_lwr(q + pp), xs_upr(q + pp)).t();
          arma::mat tempmat(n, x_grd_res);
          arma::uvec close_ind(n);
          for (arma::uword rr = 0; rr < n; rr++) {
            for (arma::uword cc = 0; cc < x_grd_res; cc++) {
              tempmat(rr, cc) = abs(prop_xs(q + pp, rr) - x_grd(q + pp, cc));
            }
          }
          for (arma::uword rr = 0; rr < n; rr++) {
            close_ind(rr) = tempmat.row(rr).index_min();
          }


          arma::uword k_x = max(z_x.row(pp)) ;
          arma::vec pi_x_ = pi_x.row(pp).t();
          prop_prior.row(q + pp) = dens_mixtnorm(prop_xs.row(q + pp).t(),
                         pi_x_(arma::span(0, k_x)),
                         mu_x(arma::span(0, k_x)),
                         sqrt(sigmasq_x(arma::span(0, k_x))),
                         xs_lwr(q + pp),
                         xs_upr(q + pp)).t();
          curr_prior.row(q + pp) = dens_mixtnorm(curr_xs.row(q + pp).t(),
                         pi_x_(arma::span(0, k_x)),
                         mu_x(arma::span(0, k_x)),
                         sqrt(sigmasq_x(arma::span(0, k_x))),
                         xs_lwr(q + pp),
                         xs_upr(q + pp)).t();
          prop_x_nz.row(q + pp) = prop_xs.row(q + pp);

          arma::mat B_basis_store_pp = B_basis_store.row(q + pp);
          prop_vars.row(q + pp) = (B_basis_store_pp.rows(close_ind) *
            exp(thetas.row(q + pp).t())).t();

          prop_pc.row(q + pp) = dtnorm(curr_xs.row(q + pp).t(),
                      prop_xs.row(q + pp).t(),
                      sd_, xs_lwr(q + pp), xs_upr(q + pp), false).t();
          prop_cp.row(q + pp) = dtnorm(prop_xs.row(q + pp).t(),
                      curr_xs.row(q + pp).t(),
                      sd_, xs_lwr(q + pp), xs_upr(q + pp), false).t();


          arma::vec prop_yxs_tmp = dist_mixtnorm(prop_xs.row(q + pp).t(),
                                                 pi_x_(arma::span(0, k_x)),
                                                 mu_x(arma::span(0, k_x)),
                                                 sqrt(
                                                   sigmasq_x(arma::span(0, k_x))
                                                   ),
                                                 xs_lwr(q + pp),
                                                 xs_upr(q + pp));
          prop_yxs_tmp.clamp(0.001, 0.999);
          for (arma::uword nn = 0; nn < n; nn++) {
            prop_yxs(q + pp, nn) = R::qnorm(prop_yxs_tmp(nn), 0, 1, true, false);
          }

          arma::vec prop_xnz_qq = prop_x_nz.row(q + pp).t();
          prop_us.row(q + pp) = ws.row(q + pp) - prop_xnz_qq(inds).t();
        }
      }

      curr_cop_prior = dmvnorm(curr_yxs.t(),
                               arma::rowvec(d, arma::fill::zeros),
                               corr_x) %
                                 exp(0.5 * sum(pow(curr_yxs, 2), 0)).t();
      prop_cop_prior = dmvnorm(prop_yxs.t(),
                               arma::rowvec(d, arma::fill::zeros),
                               corr_x) %
                                 exp(0.5 * sum(pow(prop_yxs, 2), 0)).t();

      for (arma::uword dd = 0; dd < d; dd++) {
        arma::uword k_u = max(z_u.row(dd));
        if (type_u == 0) {
          arma::uvec dd_(1, arma::fill::value(dd));
          arma::vec pi_u_dd = pi_u.row(dd).t();
          arma::vec curr_tmp = fu_mixnorm(curr_us.row(dd).t(), 0,
                                          sqrt(curr_vars.submat(dd_,inds).t()),
                                          pi_u_dd(arma::span(0, k_u)),
                                          params_u.rows(arma::span(0, k_u)));
          arma::vec prop_tmp = fu_mixnorm(prop_us.row(dd).t(), 0,
                                          sqrt(prop_vars.submat(dd_,inds).t()),
                                          pi_u_dd(arma::span(0, k_u)),
                                          params_u.rows(arma::span(0, k_u)));
          curr_marglik.row(dd) = tapply_prod(curr_tmp, mis).t();
          prop_marglik.row(dd) = tapply_prod(prop_tmp, mis).t();

          arma::vec prop_yus_tmp = dist_mixnorm(prop_us.row(dd).t(),
                                                0,
                                                sqrt(prop_vars.submat(dd_, inds).t()),
                                                pi_u_dd(arma::span(0, k_u)),
                                                params_u.rows(
                                                  arma::span(0, k_u)
                                                  ));
          prop_yus_tmp.clamp(0.001, 0.999);
          for (arma::uword NN = 0; NN < N; NN++) {
            prop_yus(dd, NN) = R::qnorm(prop_yus_tmp(NN), 0, 1, true, false);
          }
        } else if (type_u == 1) {
          arma::vec curr_tmp(N);
          arma::vec prop_tmp(N);
          arma::vec prop_yus_tmp(N);
          for (arma::uword NN = 0; NN < N; NN++) {
            curr_tmp(NN) = R::dnorm(curr_us(dd, NN), 0, sqrt(curr_vars(dd, inds(NN))), false);
            prop_tmp(NN) = R::dnorm(prop_us(dd, NN), 0, sqrt(prop_vars(dd, inds(NN))), false);
            prop_yus_tmp(NN) = R::pnorm(prop_us(dd, NN), 0, sqrt(prop_vars(dd, inds(NN))), true, false);
          }

          curr_marglik.row(dd) = tapply_prod(curr_tmp, mis).t();
          prop_marglik.row(dd) = tapply_prod(prop_tmp, mis).t();
          prop_yus_tmp.clamp(0.001, 0.999);
          for (arma::uword NN = 0; NN < N; NN++) {
            prop_yus(dd, NN) = R::qnorm(prop_yus_tmp(NN), 0, 1, true, false);
          }
        } else if (type_u == 2) {
          arma::vec curr_tmp(N);
          arma::vec prop_tmp(N);
          arma::vec prop_yus_tmp(N);
          for (arma::uword NN = 0; NN < N; NN++) {
            curr_tmp(NN) = 0.5 * R::dexp(abs(curr_us(dd, NN)), 1.0 / sqrt(curr_vars(dd, inds(NN))) / 2.0, false);
            prop_tmp(NN) = 0.5 * R::dexp(abs(prop_us(dd, NN)), 1.0 / sqrt(prop_vars(dd, inds(NN))) / 2.0, false);
            if (prop_us(dd, NN) <= 0) {
              prop_yus_tmp(NN) = 0.5 * std::exp(prop_us(dd, NN) / sqrt(prop_vars(dd, inds(NN))) / 2.0);
            } else {
              prop_yus_tmp(NN) = 1 - 0.5 * std::exp(-prop_us(dd, NN) / sqrt(prop_vars(dd, inds(NN))) / 2.0);
            }
          }

          curr_marglik.row(dd) = tapply_prod(curr_tmp, mis).t();
          prop_marglik.row(dd) = tapply_prod(prop_tmp, mis).t();
          prop_yus_tmp.clamp(0.001, 0.999);
          for (arma::uword NN = 0; NN < N; NN++) {
            prop_yus(dd, NN) = R::qnorm(prop_yus_tmp(NN), 0, 1, true, false);
          }
        }
      }

      if (q > 0) {
        for (arma::uword qq = 0; qq < q; qq++) {
          curr_marglik.row(qq) %= curr_marglik_e.row(qq);
          prop_marglik.row(qq) %= prop_marglik_e.row(qq);
        }
      }

      arma::vec curr_cop_lik_tmp = dmvnorm(curr_yus.t(),
                                           arma::rowvec(d, arma::fill::zeros),
                                           corr_e) %
                                             exp(0.5 *
                                             sum(pow(curr_yus, 2), 0)).t();
      arma::vec prop_cop_lik_tmp = dmvnorm(prop_yus.t(),
                                           arma::rowvec(d, arma::fill::zeros),
                                           corr_e) %
                                             exp(0.5 *
                                             sum(pow(prop_yus, 2), 0)).t();
      curr_cop_like = tapply_prod(curr_cop_lik_tmp, mis);
      prop_cop_like = tapply_prod(prop_cop_lik_tmp, mis);

      arma::vec mh_rat = (prod(prop_pc).t() % prod(prop_prior).t() %
        prop_cop_prior % prod(prop_marglik).t() % prop_cop_like) /
          (prod(prop_cp).t() % prod(curr_prior).t() % curr_cop_prior %
            prod(curr_marglik).t() % curr_cop_like);
      if (arma::any(arma::vectorise(mh_rat!=mh_rat))){
        mh_rat.replace(arma::datum::nan, 0);
      }

      arma::vec u(n, arma::fill::randu);
      arma::uvec inds_to_replace = arma::find(u < mh_rat);
      curr_xs.cols(inds_to_replace) = prop_xs.cols(inds_to_replace);
      xs.cols(inds_to_replace) = curr_xs.cols(inds_to_replace);
      curr_x_nz.cols(inds_to_replace) = prop_x_nz.cols(inds_to_replace); // changed
      x_nz.cols(inds_to_replace) = curr_x_nz.cols(inds_to_replace); // changed
      curr_vars.cols(inds_to_replace) = prop_vars.cols(inds_to_replace);
      vars.cols(inds_to_replace) = curr_vars.cols(inds_to_replace);
      curr_us = ws - x_nz.cols(inds);
      us = curr_us;

      if (q > 0) {
        curr_x_e.cols(inds_to_replace) = prop_x_e.cols(inds_to_replace);
        x_e.cols(inds_to_replace) = curr_x_e.cols(inds_to_replace);
        curr_u_e = w_e - x_e.cols(inds);
        u_e = curr_u_e;
      }
    }

    // Update params_u
    arma::uword k_u = max(max(z_u));
    if (it > 2000) {
      niter_u = 1;
    }
    for (arma::uword rr = 0; rr < niter_u; rr++) {
      for (arma::uword kk = 0; kk <= k_u; kk++) {
        arma::vec uspool({});
        arma::vec varspool({});
        for (arma::uword dd = 0; dd < d; dd++) {
          arma::vec u_ = us.row(dd).t();
          arma::vec var_ = vars.row(dd).t();
          arma::uvec tmp = arma::find(z_u.row(dd) == kk);
          uspool = arma::join_vert(uspool, u_(tmp));
          varspool = arma::join_vert(varspool, var_(inds(tmp)));
        }
        arma::vec prop_params_u =
          r_tnorm_proposal_params_restricted_mix_norm(params_u.row(kk).t());
        double prop_log_prior = R::dnorm(prop_params_u(1),
                                         0, sqrt(sigmasq_epsilon), true) -
                                           (a_epsilon + 1) *
                                           (log(prop_params_u(2)) +
                                           log(prop_params_u(3))) -
                                           b_epsilon *
                                           (1.0 / prop_params_u(2) +
                                           1.0 / prop_params_u(3));
        double curr_log_prior = R::dnorm(params_u(kk, 1),
                                         0, sqrt(sigmasq_epsilon), true) -
                                           (a_epsilon + 1) *
                                           (log(params_u(kk, 2)) +
                                           log(params_u(kk, 3))) -
                                           b_epsilon *
                                           (1.0 / params_u(kk, 2) +
                                           1.0 / params_u(kk, 3));
        arma::vec mu_(uspool.n_elem, arma::fill::zeros);
        arma::vec tmp_prop_loglik = log(d_restricted_mixnorm(
          uspool, mu_, sqrt(varspool), prop_params_u.t()
        ));
        arma::vec tmp_curr_loglik = log(d_restricted_mixnorm(
          uspool, mu_, sqrt(varspool), params_u.row(kk)
        ));
        double prop_loglik = sum(tmp_prop_loglik);
        double curr_loglik = sum(tmp_curr_loglik);
        double log_acc_prob = prop_log_prior + prop_loglik -
          curr_log_prior - curr_loglik;
        log_acc_prob += diff_log_curr_to_prop_tnorm_proposal_params_restricted_mix_norm(prop_params_u,params_u.row(kk).t());//changed
        if (log(R::runif(0, 1)) < log_acc_prob) {
          params_u.row(kk) = prop_params_u.t();
        }
      }
    }
    if (((k_u + 1) < K_epsilon) & (type_u == 0)) {
      for (arma::uword kk = k_u + 1; kk < K_epsilon; kk++) {
        params_u.row(kk) = r_proposal_params_restricted_mix_norm(1, 1, sigmasq_epsilon,
                     a_epsilon, b_epsilon, a_epsilon, b_epsilon).t();
      }
    }

    for (arma::uword dd = 0; dd < d; dd++) {
      // Update z_u
      if (type_u == 0) {
        for (arma::uword NN = 0; NN < N; NN++) {
          arma::vec prob_u = pi_u.row(dd).t();
          prob_u %= d_restricted_mixnorm_sclr(us(dd, NN), 0,
                                              sqrt(vars(dd, inds(NN))), params_u);
          if (arma::any(arma::vectorise(prob_u != prob_u || prob_u == arma::datum::inf || prob_u == -arma::datum::inf))){
            prob_u.replace(arma::datum::nan, 0);
            prob_u.replace(arma::datum::inf, 0);
          }
          arma::uvec val = sample(K_epsilon, 1, true, prob_u);
          z_u(dd, NN) = val(0) - 1;
        }
        n_kk_u.row(dd) = arma::hist(z_u.row(dd).t(),
                   arma::linspace<arma::uvec>(0, K_epsilon - 1, K_epsilon)).t();
        pi_u.row(dd) = rdirichlet(
          arma::conv_to<arma::vec>::from(n_kk_u.row(dd).t()) + alpha_epsilon / K_epsilon
        ).t();
      }

      // Update yus
      arma::uvec dd_(1, arma::fill::value(dd));
      arma::vec tmp_usrs(N);
      if (type_u == 0) {
        tmp_usrs = dist_mixnorm(us.row(dd).t(),
                                0, sqrt(vars.submat(dd_, inds).t()),
                                pi_u.row(dd).t(), params_u);
      } else if (type_u == 1) {
        for (arma::uword NN = 0; NN < N; NN++) {
          tmp_usrs(NN) = R::dnorm(us(dd, NN), 0, sqrt(curr_vars(dd, inds(NN))), false);
        }
      } else if (type_u == 2) {
        for (arma::uword NN = 0; NN < N; NN++) {
          tmp_usrs(NN) = R::dexp(abs(us(dd, NN)), 1.0 / sqrt(curr_vars(dd, inds(NN))) / 2.0, false);
        }
      }
      tmp_usrs.clamp(0.001, 0.999);
      for (arma::uword NN = 0; NN < N; NN++) {
        yus(dd, NN) = R::qnorm(tmp_usrs(NN), 0, 1, true, false);
        curr_yus(dd, NN) = yus(dd, NN);
      }
    }



      for (arma::uword dd = 0; dd < d; dd++) {
        // Update theta
        if (type_v) {
          arma::mat theta_cov(K_t + 1, K_t + 1);
          theta_cov.diag() = arma::vec(K_t + 1,
                        arma::fill::value(sig_tune_theta_1));
          theta_cov += sig_tune_theta_2 * prop_sig_thetas.row(dd);
          arma::vec prop_thetas_ = rmvnorm(curr_thetas.row(dd).t(),
                                          theta_cov).t();
          prop_thetas.row(dd) = prop_thetas_.t();
           prop_thetas(dd, 0) = 2 * prop_thetas(dd, 1) - prop_thetas(dd, 2); // changed

          arma::mat tempmat(n, x_grd_res);
          arma::uvec close_ind(n);
          for (arma::uword rr = 0; rr < n; rr++) {
            for (arma::uword cc = 0; cc < x_grd_res; cc++) {
              tempmat(rr, cc) = abs(x_nz(dd, rr) - x_grd(dd, cc)); // changed from prop_x_nz to x_nz
            }
          }
          for (arma::uword rr = 0; rr < n; rr++) {
            close_ind(rr) = tempmat.row(rr).index_min();
          }
          arma::mat B_basis_store_dd = B_basis_store.row(dd);
          prop_vars.row(dd) = (B_basis_store_dd.rows(close_ind) *
            exp(prop_thetas.row(dd).t())).t();

          arma::vec curr_log_prior = -curr_thetas.row(dd) * P_t *
            curr_thetas.row(dd).t() / (2.0 * sigmasq_vartheta(dd));
          arma::vec prop_log_prior = -prop_thetas.row(dd) * P_t *
            prop_thetas.row(dd).t() / (2.0 * sigmasq_vartheta(dd));

          double curr_log_lik = 0;
          double prop_log_lik = 0;
          for (arma::uword NN = 0; NN < N; NN++) {
            arma::vec cval(1);
            arma::vec pval(1);
            if (type_u == 0) {
              cval = d_restricted_mixnorm_sclr(us(dd, NN), 0,
                                              sqrt(
                                                curr_vars(dd, inds(NN))
                                              ),
                                              params_u.row(z_u(dd, NN)));
              pval = d_restricted_mixnorm_sclr(us(dd, NN), 0,
                                              sqrt(
                                                prop_vars(dd, inds(NN))
                                              ),
                                              params_u.row(z_u(dd, NN)));
            } else if (type_u == 1) {
              cval(0) = R::dnorm(us(dd, NN), 0, sqrt(curr_vars(dd, inds(NN))), false);
              pval(0) = R::dnorm(us(dd, NN), 0, sqrt(prop_vars(dd, inds(NN))), false);
            } else if (type_u == 2) {
              cval(0) = 0.5 * R::dexp(abs(us(dd, NN)), 1.0 / sqrt(curr_vars(dd, inds(NN))) / 2.0, false);
              pval(0) = 0.5 * R::dexp(abs(us(dd, NN)), 1.0 / sqrt(prop_vars(dd, inds(NN))) / 2.0, false);
            }
            if (cval(0) == 0) {
              cval(0) = 0.0000000001;
            }
            if (pval(0) == 0) {
              pval(0) = 0.0000000001;
            }
            curr_log_lik += log(cval(0));
            prop_log_lik += log(pval(0));
          }

          double log_mh_rat = prop_log_prior(0) + prop_log_lik -
            curr_log_prior(0) - curr_log_lik;
          if (log_mh_rat == arma::datum::nan) {
            log_mh_rat = -arma::datum::inf;
          }
          if (log(R::runif(0, 1)) < log_mh_rat) {
            curr_thetas.row(dd) = prop_thetas.row(dd);
            thetas.row(dd) = curr_thetas.row(dd);
            curr_vars.row(dd) = prop_vars.row(dd);
            vars.row(dd) = curr_vars.row(dd);
          }
        }
        // Update sigmasq_vartheta
        arma::vec sigmasq_vartheta_rate = b_vartheta + thetas.row(dd) * P_t * thetas.row(dd).t();
        sigmasq_vartheta(dd) = 1.0 / R::rgamma(a_vartheta + (K_t + 1) / 2.0, 1.0 / sigmasq_vartheta_rate(0));
      }


    if (p > 0) {
      for (arma::uword pp = 0; pp < p; pp++) {
        arma::vec xs_trans_tmp_num(n);
        arma::vec xs_trans_tmp_den(n);
        for (arma::uword nn = 0; nn < n; nn++) {
          arma::uword zz = z_x(pp, nn);
          double mu_ = mu_x(zz);
          double sigsq_ = sigmasq_x(zz);
          xs_trans_tmp_num(nn) = R::pnorm((xs(q + pp, nn) - mu_) / sqrt(sigsq_),
                           0, 1, true, false);
          xs_trans_tmp_num(nn) -= R::pnorm((xs_lwr(q + pp) - mu_) /
            sqrt(sigsq_), 0, 1, true, false);
          xs_trans_tmp_den(nn) = R::pnorm((xs_upr(q + pp) - mu_) / sqrt(sigsq_),
                           0, 1, true, false);
          xs_trans_tmp_den(nn) -= R::pnorm((xs_lwr(q + pp) - mu_) /
            sqrt(sigsq_), 0, 1, true, false);
          double xs_trans_tmp = xs_trans_tmp_num(nn) / xs_trans_tmp_den(nn);
          if (xs_trans_tmp >= 1) {
            xs_trans_tmp = 0.999;
          }
          if (xs_trans_tmp <= 0) {
            xs_trans_tmp = 0.001;
          }
          xs_trans_tmp = mu_ + sqrt(sigsq_) * R::qnorm(xs_trans_tmp, 0, 1, true,
                                    false);
          xs_trans(pp, nn) = xs_trans_tmp;
        }
      }
    }

    if (q > 0) {
      for (arma::uword qq = 0; qq < q; qq++) {
        // Update w_e
        arma::vec x_e_ = x_e.row(qq).t();
        arma::vec xs_star = x_e_(inds);
        arma::uvec inc_tmp = arma::find(inc.row(qq) == 1);
        for (arma::uword jj = 0; jj < inc_tmp.n_elem; jj++) {
          w_e(qq, inc_tmp(jj)) = rtnorm(xs_star(inc_tmp(jj)), 1,
              -arma::datum::inf, 0);
        }
        arma::uvec ic_tmp = arma::find(ic.row(qq) == 1);
        for (arma::uword jj = 0; jj < ic_tmp.n_elem; jj++) {
          w_e(qq, ic_tmp(jj)) = rtnorm(xs_star(ic_tmp(jj)), 1,
              0, arma::datum::inf);
        }
        // Update theta_e
        arma::mat tempmat(n, x_grd_res);
        arma::uvec close_ind(n);
        for (arma::uword rr = 0; rr < n; rr++) {
          for (arma::uword cc = 0; cc < x_grd_res; cc++) {
            tempmat(rr, cc) = abs(xs(qq, rr) - x_grd(qq, cc));
          }
        }
        for (arma::uword rr = 0; rr < n; rr++) {
          close_ind(rr) = tempmat.row(rr).index_min();
        }
        arma::vec mu_theta_e(K_t + 1, arma::fill::zeros);
        arma::mat sigma_theta_e(K_t + 1, K_t + 1, arma::fill::zeros);
        arma::mat B_basis_store_qq = B_basis_store.row(qq);
        for (arma::uword nn = 0; nn < n; nn++) {
          for (arma::uword hh = inds1(nn); hh <= inds2(nn); hh++) {
            mu_theta_e += w_e(qq, hh - 1) *
              B_basis_store_qq.row(close_ind(nn)).t();
          }
          sigma_theta_e += mis(nn) * B_basis_store_qq.row(close_ind(nn)).t() *
            B_basis_store_qq.row(close_ind(nn));
        }
        sigma_theta_e = arma::inv(P_t / sigmasq_beta(qq) +
          arma::inv(sigma0_thetas_eps) + sigma_theta_e);
        mu_theta_e = sigma_theta_e * (arma::inv(sigma0_thetas_eps) *
          mu0_thetas_eps.row(qq).t() + mu_theta_e);
        arma::vec cond_mu = cond * mu_theta_e;
        arma::vec proposed_theta_e = mu_theta_e +
          icond * rmvtnorm(cond_lwr - cond_mu,
                           cond_upr - cond_mu,
                           cond * sigma_theta_e * cond.t()); // changed

        prop_x_e.row(qq) = (B_basis_store_qq.rows(close_ind) *
          proposed_theta_e).t();
        for (arma::uword nn = 0; nn < n; nn++) {
          prop_x_nz(qq, nn) = xs(qq, nn) / R::pnorm(prop_x_e(qq, nn),
                    0, 1, true, false);
        }
        arma::vec prop_xnz_qq = prop_x_nz.row(qq).t();
        prop_us.row(qq) = ws.row(qq) - prop_xnz_qq(inds).t();
        double curr_log_lik = 0;
        double prop_log_lik = 0;
        for (arma::uword NN = 0; NN < N; NN++) {
          arma::vec cval(1);
          arma::vec pval(1);
          if (type_u == 0) {
            cval = d_restricted_mixnorm_sclr(prop_us(qq, NN), 0,
                                             sqrt(
                                               vars(qq, inds(NN))
                                             ),
                                             params_u.row(z_u(qq, NN)));
            pval = d_restricted_mixnorm_sclr(curr_us(qq, NN), 0,
                                             sqrt(
                                               vars(qq, inds(NN))
                                             ),
                                             params_u.row(z_u(qq, NN)));

          } else if (type_u == 1) {
            cval(0) = R::dnorm(prop_us(qq, NN), 0, sqrt(vars(qq, inds(NN))), false);
            pval(0) = R::dnorm(curr_us(qq, NN), 0, sqrt(vars(qq, inds(NN))), false);
          } else if (type_u == 2) {
            cval(0) = 0.5 * R::dexp(abs(prop_us(qq, NN)), 1.0 / sqrt(vars(qq, inds(NN))) / 2.0, false);
            pval(0) = 0.5 * R::dexp(abs(curr_us(qq, NN)), 1.0 / sqrt(vars(qq, inds(NN))) / 2.0, false);
          }
          if (cval(0) == 0) {
            cval(0) = 0.0000000001;
          }
          if (pval(0) == 0) {
            pval(0) = 0.0000000001;
          }
          curr_log_lik += log(cval(0));
          prop_log_lik += log(pval(0));
        }
        double log_mh_rat = prop_log_lik - curr_log_lik;
        if (log_mh_rat == arma::datum::nan) {
          log_mh_rat = -arma::datum::inf;
        }
        if (log(R::runif(0, 1)) < log_mh_rat) {
          theta_e.row(qq) = proposed_theta_e.t();
        }

        curr_x_e.row(qq) = (B_basis_store_qq.rows(close_ind) *
          theta_e.row(qq).t()).t();
        x_e.row(qq) = curr_x_e.row(qq);
        for (arma::uword nn = 0; nn < n; nn++) {
          curr_x_nz(qq, nn) = xs(qq, nn) / R::pnorm(x_e(qq, nn),
                    0, 1, true, false);
          x_nz(qq, nn) = curr_x_nz(qq, nn);
        }
        x_nz.row(qq).clamp(min(knots_t.row(qq)), max(knots_t.row(qq)) - 0.001);

        for (arma::uword rr = 0; rr < n; rr++) {
          for (arma::uword cc = 0; cc < x_grd_res; cc++) {
            tempmat(rr, cc) = abs(x_nz(qq, rr) - x_grd(qq, cc));
          }
        }
        for (arma::uword rr = 0; rr < n; rr++) {
          close_ind(rr) = tempmat.row(rr).index_min();
        }
        curr_vars.row(qq) = (B_basis_store_qq.rows(close_ind) *
          exp(thetas.row(qq).t())).t();
        vars.row(qq) = curr_vars.row(qq);
        // Update s2t_e
        arma::vec sigmasq_beta_rate = b_beta(qq) + theta_e.row(qq) * P_t *
          theta_e.row(qq).t();
        sigmasq_beta(qq) = 1.0 / R::rgamma(a_beta(qq) + (K_t + 1) / 2.0,
                1.0 / sigmasq_beta_rate(0));
        // Update ws when not observed
        arma::vec mu_w(N);
        arma::vec sigmasq_w(N);
        if (type_u == 0) {
          arma::vec c2_u(N);
          for (arma::uword NN = 0; NN < N; NN++) {
            c2_u(NN) = R::rbinom(1, params_u(z_u(qq, NN)));
          }
          arma::mat mu_u(params_u.n_rows, 2);
          mu_u.col(0) = params_u.col(1) % (1.0 - params_u.col(0)) /
            (params_u.col(0) % params_u.col(0) +
              (1.0 - params_u.col(0)) % (1.0 - params_u.col(0)));
          mu_u.col(1) = -mu_u.col(0) % params_u.col(0) / (1.0 - params_u.col(0));
          for (arma::uword NN = 0; NN < N; NN++) {
            mu_w(NN) = x_nz(qq, inds(NN));
            mu_w(NN) += sqrt(vars(qq, inds(NN))) * mu_u(z_u(qq, NN), c2_u(NN));
            sigmasq_w(NN) = vars(qq, inds(NN)) * params_u(z_u(qq, NN),
                      c2_u(NN) + 2);
          }
        } else {
          for (arma::uword NN = 0; NN < N; NN++) {
            mu_w(NN) = x_nz(qq, inds(NN));
            sigmasq_w(NN) = vars(qq, inds(NN));
          }
        }

        for (arma::uword idx = 0; idx < inc_tmp.n_elem; idx++) {
          ws(qq, inc_tmp(idx)) = R::rnorm(mu_w(inc_tmp(idx)),
             sqrt(sigmasq_w(inc_tmp(idx))));
        }
      }
      // Re-update u_e and us
      curr_u_e = w_e - x_e.cols(inds);
      u_e = curr_u_e;
      curr_us = ws - x_nz.cols(inds);
      us = curr_us;
    }

    // Rcpp::Rcout << "\r\tRecord output";
    if ((it >= nburn) & ((it - nburn) % nthin == 0)) {
      arma::uword idx = (it - nburn) / nthin;
      ws_mcmc.row(idx) = ws;
      xs_mcmc.row(idx) = xs;
      z_x_mcmc.row(idx) = z_x;
      pi_x_mcmc.row(idx) = pi_x;
      for (arma::uword pp = 0; pp < p; pp++){
        for (arma::uword kk = 0; kk < K_x; kk++){
          mu_x_mcmc(idx,pp,kk) = mu_x(kk);
          sig_x_mcmc(idx,pp,kk) = sqrt(sigmasq_x(kk));
        }
      }
      x_e_mcmc.row(idx) = x_e;
      x_nz_mcmc.row(idx) = x_nz;
      thetas_x_dens_mcmc.row(idx ) = thetas_x_density;
      thetas_mcmc.row(idx) = thetas;
      theta_e_mcmc.row(idx) = theta_e;


      if (type_u == 0) {
        for (arma::uword dd = 0; dd < d; dd++) {
          arma::uword k_u = max(z_u.row(dd));
          arma::vec pi_u_ = pi_u.row(dd).t();
          var_es(dd) = var_e_fn(pi_u_(arma::span(0, k_u)),
                            params_u.rows(arma::span(0, k_u)));
        }
      } else {
        var_es = arma::vec(d, arma::fill::ones);
      }
      if (!type_v) {
        for (arma::uword dd = 0; dd < d; dd++) {
          var_es(dd) *= 1.0 / R::rgamma(1 + 0.5 * N, 1.0 / (1 + 0.5 * arma::accu(us.row(dd) % us.row(dd))));
        }
      }
      var_es_mcmc.row(idx) = var_es.t();

      us_mcmc.row(idx) = us;
      z_u_mcmc.row(idx) = z_u;
      k_u_mcmc.row(idx) = max(z_u, 1).t();
      pi_u_mcmc.row(idx) = pi_u;
      p_u_mcmc.row(idx) = params_u.col(0).t();
      mu_u_mcmc.row(idx) = params_u.col(1).t();
      sig1_u_mcmc.row(idx) = sqrt(params_u.col(2).t());
      sig2_u_mcmc.row(idx) = sqrt(params_u.col(3).t());

      corr_x_mcmc.row(idx) = corr_x;
      corr_e_mcmc.row(idx) = corr_e;
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
  // Sublist to avoid 20 element cap for Rcpp::List obj
  Rcpp::List u_res = Rcpp::List::create(Rcpp::_("u") = us_mcmc,
                                        Rcpp::_("z_u") = z_u_mcmc + 1,
                                        Rcpp::_("k_u") = k_u_mcmc + 1,
                                        Rcpp::_("pi_u") = pi_u_mcmc,
                                        Rcpp::_("p_u") = p_u_mcmc,
                                        Rcpp::_("mu_u") = mu_u_mcmc,
                                        Rcpp::_("sig1_u") = sig1_u_mcmc,
                                        Rcpp::_("sig2_u") = sig2_u_mcmc);
  Rcpp::List corr_res = Rcpp::List::create(Rcpp::_("corr_x") = corr_x_mcmc,
                                           Rcpp::_("corr_e") = corr_e_mcmc);

  return Rcpp::List::create(Rcpp::_("w") = ws_mcmc,
                            Rcpp::_("x") = xs_mcmc,
                            Rcpp::_("z_x") = z_x_mcmc + 1,
                            Rcpp::_("pi_x") = pi_x_mcmc,
                            Rcpp::_("mu_x") = mu_x_mcmc,
                            Rcpp::_("sig_x") = sig_x_mcmc,
                            Rcpp::_("theta_x") = thetas_x_dens_mcmc,
                            Rcpp::_("theta_v") = thetas_mcmc,
                            Rcpp::_("x_nz") = x_nz_mcmc,
                            Rcpp::_("x_e") = x_e_mcmc,
                            Rcpp::_("theta_e") = theta_e_mcmc,
                            Rcpp::_("var_e") = var_es_mcmc,
                            Rcpp::_("u_res") = u_res,
                            Rcpp::_("corr_res") = corr_res);
}
