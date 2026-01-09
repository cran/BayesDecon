#include <RcppArmadillo.h>

// Lite port of truncated normal functions from `msm`, Jackson (2021)

// TODO: Parallelize, reduce calls to pnorm upr/lwr
// [[Rcpp::export]]
double dtnorm_sclr(double x, double mean, double sd, double lwr, double upr,
                   bool log) {
  if (upr < lwr || sd < 0) {
    Rcpp::stop("dtnorm: invalid parameters");
  }
  if (x < lwr || x > upr) {
    if (log) {
      return -arma::datum::inf;
    } else {
      return 0;
    }
  } else {
    double denom = R::pnorm(upr, mean, sd, true, false) -
      R::pnorm(lwr, mean, sd, true, false);
    double xtmp = R::dnorm(x, mean, sd, log);
    if (log) {
      xtmp -= std::log(denom);
    } else {
      xtmp /= denom;
    }
    return xtmp;
  }
}

// [[Rcpp::export]]
arma::vec dtnorm(arma::vec x, double mean, double sd, double lwr, double upr,
                 bool log) {
  if (upr < lwr || sd < 0) {
    Rcpp::stop("dtnorm: invalid parameters");
  }
  arma::vec ret(x.n_elem, arma::fill::zeros);
  arma::uvec idx_lwr = arma::find(x < lwr);
  arma::uvec idx_upr = arma::find(x > upr);
  if (log) {
    ret(idx_lwr) = arma::vec(idx_lwr.n_elem,
        arma::fill::value(-arma::datum::inf));
    ret(idx_upr) = arma::vec(idx_upr.n_elem,
        arma::fill::value(-arma::datum::inf));
  } else {
    ret(idx_lwr) = arma::vec(idx_lwr.n_elem, arma::fill::zeros);
    ret(idx_upr) = arma::vec(idx_upr.n_elem, arma::fill::zeros);
  }
  arma::uvec ind = arma::find((x >= lwr) && (x <= upr));
  if (ind.n_elem > 0) {
    double denom = R::pnorm(upr, mean, sd, true, false) -
      R::pnorm(lwr, mean, sd, true, false);
    arma::vec xtmp(x.n_elem);
    for (arma::uword i = 0; i < x.n_elem; i++) {
      xtmp(i) = R::dnorm(x(i), mean, sd, log);
    }
    if (log) {
      xtmp -= std::log(denom);
    } else {
      xtmp /= denom;
    }
    ret(ind) = xtmp(ind);
  }
  return ret;
}

arma::vec dtnorm(arma::vec x, arma::vec mean, double sd, double lwr, double upr,
                 bool log) {
  arma::vec ret(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; i++) {
    ret(i) = dtnorm_sclr(x(i), mean(i), sd, lwr, upr, log);
  }
  return ret;
}

arma::vec dtnorm(arma::vec x, arma::vec mean, arma::vec sd, double lwr, double upr,
                 bool log) {
  arma::vec ret(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; i++) {
    ret(i) = dtnorm_sclr(x(i), mean(i), sd(i), lwr, upr, log);
  }
  return ret;
}

// [[Rcpp::export]]
arma::vec ptnorm(arma::vec q, double mean, double sd, double lwr, double upr,
                 bool lower_tail, bool log_p) {
  if (upr < lwr) {
    Rcpp::stop("ptnorm: upr < lwr");
  }
  arma::vec ret(q.n_elem);
  arma::uvec idx_lwr = arma::find(q < lwr);
  arma::uvec idx_upr = arma::find(q > upr);
  if (lower_tail) {
    ret(idx_lwr) = arma::vec(idx_lwr.n_elem, arma::fill::zeros);
    ret(idx_upr) = arma::vec(idx_upr.n_elem, arma::fill::ones);
  } else {
    ret(idx_lwr) = arma::vec(idx_lwr.n_elem, arma::fill::ones);
    ret(idx_upr) = arma::vec(idx_upr.n_elem, arma::fill::zeros);
  }
  arma::uvec ind = arma::find((q >= lwr) && (q <= upr));
  if (ind.n_elem > 0) {
    double denom = R::pnorm(upr, mean, sd, true, false) -
      R::pnorm(lwr, mean, sd, true, false);
    arma::vec qtmp(q.n_elem);
    for (arma::uword i = 0; i < q.n_elem; i++) {
      if (lower_tail) {
        qtmp(i) = R::pnorm(q(i), mean, sd, true, false) -
          R::pnorm(lwr, mean, sd, true, false);
      } else {
        qtmp(i) = R::pnorm(upr, mean, sd, true, false) -
          R::pnorm(q(i), mean, sd, true, false);
      }
    }
    if (log_p) {
      qtmp = arma::log(qtmp) - std::log(denom);
    } else {
      qtmp /= denom;
    }
    ret(ind) = qtmp(ind);
  }
  return ret;
}

arma::vec ptnorm(arma::vec q, arma::vec mean, arma::vec sd, double lwr, double upr,
                 bool lower_tail, bool log_p) {
  if (upr < lwr) {
    Rcpp::stop("ptnorm: upr < lwr");
  }
  arma::vec ret(q.n_elem);
  arma::uvec idx_lwr = arma::find(q < lwr);
  arma::uvec idx_upr = arma::find(q > upr);
  if (lower_tail) {
    ret(idx_lwr) = arma::vec(idx_lwr.n_elem, arma::fill::zeros);
    ret(idx_upr) = arma::vec(idx_upr.n_elem, arma::fill::ones);
  } else {
    ret(idx_lwr) = arma::vec(idx_lwr.n_elem, arma::fill::ones);
    ret(idx_upr) = arma::vec(idx_upr.n_elem, arma::fill::zeros);
  }
  arma::uvec ind = arma::find((q >= lwr) && (q <= upr));
  if (ind.n_elem > 0) {
    arma::vec denom(q.n_elem);
    for (arma::uword i = 0; i < q.n_elem; i++) {
      denom(i) = R::pnorm(upr, mean(i), sd(i), true, false) -
        R::pnorm(lwr, mean(i), sd(i), true, false);
    }
    arma::vec qtmp(q.n_elem);
    for (arma::uword i = 0; i < q.n_elem; i++) {
      if (lower_tail) {
        qtmp(i) = R::pnorm(q(i), mean(i), sd(i), true, false) -
          R::pnorm(lwr, mean(i), sd(i), true, false);
      } else {
        qtmp(i) = R::pnorm(upr, mean(i), sd(i), true, false) -
          R::pnorm(q(i), mean(i), sd(i), true, false);
      }
    }
    if (log_p) {
      qtmp = arma::log(qtmp) - arma::log(denom);
    } else {
      qtmp /= denom;
    }
    ret(ind) = qtmp(ind);
  }
  return ret;
}

// Robert (Stat. Comp (1995), 5, 121-5), cited by Jackson (2021)
// [[Rcpp::export]]
arma::vec rtnorm(arma::uword n, double mean, double sd, double lwr,
                 double upr) {
  if (upr < lwr) {
    Rcpp::stop("ptnorm: upr < lwr");
  }
  lwr = (lwr - mean) / sd;
  upr = (upr - mean) / sd;
  arma::uvec ind = arma::linspace<arma::uvec>(1, n, n);
  arma::vec ret(n);
  if ((lwr < 0 && upr == arma::datum::inf) ||
      (lwr == -arma::datum::inf && upr > 0) ||
      (lwr < 0 && upr > 0 && upr - lwr > sqrt(2 * 3.14159))) {
    for (arma::uword i = 0; i < n; i++) {
      double y = R::rnorm(0, 1);
      while ((y < lwr) | (y > upr)) {
        y = R::rnorm(0, 1);
      }
      ret(i) = y;
    }
    return ret * sd + mean;
  } else if ((lwr >= 0) & (upr > lwr + 2.0 * sqrt(exp(1.0)) /
    (lwr + sqrt(lwr * lwr + 4)) * exp((lwr * lwr - lwr * sqrt(lwr * lwr + 4)) /
    4))) {
    for (arma::uword i = 0; i < n; i++) {
      double a = (lwr + sqrt(lwr * lwr + 4)) / 2;
      double z = R::rexp(1 / a) + lwr;
      double u = R::runif(0, 1);
      while (!((u <= exp(-(z - a) * (z - a) / 2)) && (z <= upr))) {
        z = R::rexp(1 / a) + lwr;
        u = R::runif(0, 1);
      }
      ret(i) = z;
    }
    return ret * sd + mean;
  } else if ((upr <= 0) & (-lwr > -upr + 2.0 * sqrt(exp(1.0)) /
    (-upr + sqrt(upr * upr + 4)) * exp((upr * 2 - -upr * sqrt(upr * upr + 4)) /
    4))) {
    for (arma::uword i = 0; i < n; i++) {
      double a = (-upr + sqrt(upr * upr + 4)) / 2;
      double z = R::rexp(1 / a) - upr;
      double u = R::runif(0, 1);
      while (!((u <= exp(-(z - a) * (z - a) / 2)) && (z <= -lwr))) {
        z = R::rexp(1 / a) - upr;
        u = R::runif(0, 1);
      }
      ret(i) = -z;
    }
    return ret * sd + mean;
  } else {
    for (arma::uword i = 0; i < n; i++) {
      double z = R::runif(lwr, upr);
      double rho;
      if (lwr < 0) {
        rho = exp((lwr * lwr - z * z) / 2);
      } else if (upr < 0) {
        rho = exp((upr * upr - z * z) / 2);
      } else {
        rho = exp(-z * z / 2);
      }
      double u = R::runif(0, 1);
      while (u > rho) {
        z = R::runif(lwr, upr);
        if (lwr < 0) {
          rho = exp((lwr * lwr - z * z) / 2);
        } else if (upr < 0) {
          rho = exp((upr * upr - z * z) / 2);
        } else {
          rho = exp(-z * z / 2);
        }
        u = R::runif(0, 1);
      }
      ret(i) = z;
    }
    return ret * sd + mean;
  }
}

double rtnorm(double mean, double sd, double lwr, double upr) {
  arma::vec ret = rtnorm(1, mean, sd, lwr, upr);
  return ret(0);
}

arma::vec rtnorm(arma::vec mean, double sd, double lwr, double upr) {
  arma::vec ret(mean.n_elem);
  for (arma::uword i = 0; i < mean.n_elem; i++) {
    ret(i) = rtnorm(mean(i), sd, lwr, upr);
  }
  return ret;
}

arma::vec dist_mixtnorm(arma::vec x, arma::vec pi, arma::vec mean, arma::vec sd,
                        double lwr, double upr) {
  arma::mat y(pi.n_elem, x.n_elem, arma::fill::zeros);
  for (arma::uword kk = 0; kk < pi.n_elem; kk++) {
    y.row(kk) = pi(kk) * ptnorm(x, mean(kk), sd(kk), lwr, upr, true, false).t();
  }
  return(sum(y, 0).t());
}

arma::vec dens_mixtnorm(arma::vec x, arma::vec pi, arma::vec mean, arma::vec sd,
                        double lwr, double upr) {
  arma::mat y(pi.n_elem, x.n_elem, arma::fill::zeros);
  for (arma::uword kk = 0; kk < pi.n_elem; kk++) {
    y.row(kk) = pi(kk) * dtnorm(x, mean(kk), sd(kk), lwr, upr, false).t();
  }
  return(sum(y, 0).t());
}

// Method of Kotecha & Djuric 1999 and Breslaw 1994 restricted to (0, 1)^n
// [[Rcpp::export]]
arma::vec rmvtnorm(arma::vec alpha, arma::vec beta, arma::mat cov) {
  arma::uword n = alpha.n_elem;
  arma::vec mean(alpha.n_elem, arma::fill::zeros);
  arma::vec ret(n);
  arma::uvec idxs = arma::linspace<arma::uvec>(0, n - 1, n);
  arma::mat muset(n - 1, n);
  arma::cube cov12set(n, 1, n - 1);
  arma::cube cov22invset(n, n - 1, n - 1);
  for (arma::uword i = 0; i < n; i++) {
    ret(i) = rtnorm(mean(i), sqrt(cov(i, i)), alpha(i), beta(i));
    arma::uvec idxs_ = idxs;
    idxs_.shed_row(i);
    muset.col(i) = mean(idxs_);
    cov12set.row(i) = cov.submat(arma::uvec({i}), idxs_);
    cov22invset.row(i) = arma::inv(cov.submat(idxs_, idxs_));
  }
  arma::vec mu(n - 1);
  arma::rowvec cov12(n - 1);
  arma::mat cov22inv(n - 1, n - 1);
  for (arma::uword j = 0; j < n * 5; j++) {
    for (arma::uword i = 0; i < n; i++) {
      mu = muset.col(i);
      cov12 = cov12set.row(i);
      cov22inv = cov22invset.row(i);
      arma::vec ret_ = ret;
      ret_.shed_row(i);
      arma::mat mean_cond_ = cov12 * cov22inv * (ret_ - mu);
      double mean_cond = mean(i) + mean_cond_(0, 0);
      arma::mat sd_cond_ = cov(i, i) - cov12 * cov22inv * cov12.t();
      double sd_cond = sqrt(sd_cond_(0, 0));
      ret(i) = rtnorm(mean_cond, sd_cond, alpha(i), beta(i));
    }
  }
  return ret;
}

// Monotonically increasing version of rmvtnorm
// [[Rcpp::export]]
arma::vec rmvtnorm_mono(arma::vec mean, arma::mat cov) {
  arma::uword n = mean.n_elem;
  arma::vec ret(n, arma::fill::value(arma::datum::nan));
  arma::uvec idxs = arma::linspace<arma::uvec>(0, n - 1, n);
  arma::mat muset(n - 1, n);
  arma::cube cov12set(n, 1, n - 1);
  arma::cube cov22invset(n, n - 1, n - 1);
  double lwr;
  double upr;
  for (arma::uword i = 0; i < n; i++) {
    lwr = -arma::datum::inf;
    upr = arma::datum::inf;
    if (idxs(i) > 0) {
      lwr = ret(arma::span(0, idxs(i) - 1)).max();
    }
    if (idxs(i) < (n - 1)) {
      upr = ret(arma::span(idxs(i) + 1, n - 1)).min();
    }
    ret(idxs(i)) = rtnorm(mean(idxs(i)), sqrt(cov(idxs(i), idxs(i))),
                          lwr, upr);
    arma::uvec idxs_ = idxs;
    idxs_.shed_row(i);
    muset.col(idxs(i)) = mean(idxs_);
    cov12set.row(idxs(i)) = cov.submat(arma::uvec({idxs(i)}), idxs_);
    cov22invset.row(idxs(i)) = arma::inv(cov.submat(idxs_, idxs_));
  }
  arma::vec mu(n - 1);
  arma::rowvec cov12(n - 1);
  arma::mat cov22inv(n - 1, n - 1);
  for (arma::uword j = 0; j < n * 5; j++) {
    for (arma::uword i = 0; i < n; i++) {
      mu = muset.col(idxs(i));
      cov12 = cov12set.row(idxs(i));
      cov22inv = cov22invset.row(idxs(i));
      arma::vec ret_ = ret;
      ret_.shed_row(idxs(i));
      arma::mat mean_cond_ = cov12 * cov22inv * (ret_ - mu);
      double mean_cond = mean(idxs(i)) + mean_cond_(0, 0);
      arma::mat sd_cond_ = cov(idxs(i), idxs(i)) - cov12 * cov22inv * cov12.t();
      double sd_cond = sqrt(sd_cond_(0, 0));
      lwr = -arma::datum::inf;
      upr = arma::datum::inf;
      if (idxs(i) > 0) {
        lwr = ret(arma::span(0, idxs(i) - 1)).max();
      }
      if (idxs(i) < (n - 1)) {
        upr = ret(arma::span(idxs(i) + 1, n - 1)).min();
      }
      ret(idxs(i)) = rtnorm(mean_cond, sd_cond, lwr, upr);
    }
  }
  return ret;
}

