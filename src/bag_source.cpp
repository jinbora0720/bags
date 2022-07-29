#include "matrix_distance.h"
#include "bag_source.h"
#include <RcppArmadilloExtensions/sample.h>

////////////////
// create_RnH //
////////////////
// [[Rcpp::export]]
Rcpp::List createRnH(const arma::mat& coords,
                     const Rcpp::List& idx,
                     const Rcpp::List& pptts_list,
                     const Rcpp::CharacterVector& ptts_proto,
                     const double& a, const double& c, const double& kappa,
                     bool produce_R, const Rcpp::CharacterVector& directions) {

  // coords: coordinates with partition
  // idx: index of all partitions
  // pptts_list: list of parent partitions and their idx for each partition
  // ptts_proto: prototypical partitions e.g. (*, !, 1), (*, !, 2)
  // a: temporal decay
  // c: spatial decay
  // kappa: interaction between space and time
  // produce_R = FALSE: do not save R

  std::string m;
  std::string h;
  arma::uvec ColInd = {0, 1, 2};
  Rcpp::List RnH_list;

  for (int i=0; i < ptts_proto.length(); i++) {
    // level 1: prototype partition
    m = ptts_proto(i);
    arma::uvec idx_tmp = Rcpp::as<arma::uvec>(idx[m]);
    int n_ij = idx_tmp.n_elem;
    arma::uvec onevec = arma::ones<arma::uvec>(n_ij);
    arma::uvec idx_ij = idx_tmp - onevec;
    arma::mat coords_m = coords.submat(idx_ij, ColInd);
    Rcpp::List pptts_list_byptt = pptts_list[m];
    Rcpp::List RnH_wind;

    for (int j=0; j < directions.length(); j++) {
      // level 2: wind direction
      h = directions(j);
      Rcpp::Nullable<Rcpp::List> pptts_list_byptt_bydir_ = pptts_list_byptt[h];
      Rcpp::List RnH;

      if (pptts_list_byptt_bydir_.isNotNull()) {
        // at least one parent
        Rcpp::List pptts_list_byptt_bydir(pptts_list_byptt_bydir_);
        arma::uvec idx_tmpp = Rcpp::as<arma::uvec>(pptts_list_byptt_bydir["pidx"]);
        int n_pij = idx_tmpp.n_elem;
        arma::uvec onevec2 = arma::ones<arma::uvec>(n_pij);
        arma::uvec idx_pij = idx_tmpp - onevec2;
        arma::mat coords_mh = join_cols(coords_m, coords.submat(idx_pij, ColInd));

        arma::mat spatdist = distmat_Euclidean(coords_mh.submat(0, 0, n_ij + n_pij -1, 1));
        arma::mat timedist = distmat_abs(coords_mh.submat(0, 2, n_ij + n_pij -1, 2));
        arma::mat aup1 = a*timedist+1;
        arma::mat Cor = exp(-c*spatdist/pow(aup1, kappa/2))/aup1;
        arma::mat Cipi = Cor.submat(0, n_ij, n_ij-1, n_ij+n_pij-1);
        arma::mat CinvC = Cipi*
          inv_sympd(Cor.submat(n_ij, n_ij, n_ij+n_pij-1, n_ij+n_pij-1));
        arma::mat R = Cor.submat(0, 0, n_ij-1, n_ij-1) - CinvC*Cipi.t();
        R = symmatu(R);                                                         // had to force symmetry.

        RnH["H"] = CinvC;
        if (R.is_sympd()) {
          RnH["tB"] = arma::inv(arma::chol(R, "lower"));
        } else {
          RnH["tB"] = arma::inv(arma::chol(R + 0.005*arma::eye(n_ij, n_ij), "lower"));
        }
        if (produce_R) {
          RnH["R"] = R;
        }
      } else {
        // no parents
        arma::mat spatdist = distmat_Euclidean(coords_m.submat(0, 0, n_ij -1, 1));
        arma::mat timedist = distmat_abs(coords_m.submat(0, 2, n_ij -1, 2));
        arma::mat aup1 = a*timedist+1;
        arma::mat Cor = exp(-c*spatdist/pow(aup1, kappa/2))/aup1;

        if (Cor.is_sympd()) {
          RnH["tB"] = arma::inv(arma::chol(Cor, "lower"));
        } else {
          RnH["tB"] = arma::inv(arma::chol(Cor + 0.005*arma::eye(n_ij, n_ij), "lower"));
        }
        if (produce_R) {
          RnH["R"] = Cor;
        }
      }

      RnH_wind[h] = RnH;
    }

    RnH_list[m] = RnH_wind;
  }

  return RnH_list;
}

/////////////////
// wsig_update //
/////////////////
// [[Rcpp::export]]
Rcpp::List wsigupdate(const arma::vec& res, arma::vec& w, double& sig_sq,
                      const double& tau_sq, const Rcpp::CharacterVector& z,
                      const Rcpp::List& RnH_list, const Rcpp::CharacterVector& ptts,
                      const Rcpp::List& idx, const Rcpp::List& pptts_list,
                      const Rcpp::CharacterVector& ptts_proto, const double& nd,
                      const double& as, const double& bs) {

  std::string m;
  std::string h;
  int M = ptts.length();
  Rcpp::CharacterMatrix ppartition(M, 2);

  // attach ppartition to ptt and z
  for (int i=0; i < M; i++) {
    m = ptts(i);
    h = z(i);
    Rcpp::List pptts_list_byptt = pptts_list[m];
    Rcpp::Nullable<Rcpp::List> pptts_list_byptt_bydir_ = pptts_list_byptt[h];
    if (pptts_list_byptt_bydir_.isNotNull()) {
      // at least one parent
      Rcpp::List pptts_list_byptt_bydir(pptts_list_byptt_bydir_);
      Rcpp::CharacterVector pptttmp = pptts_list_byptt_bydir["ppartition"];
      ppartition(i,Rcpp::_) = pptttmp;
    } else {
      // no parents
      ppartition(i,Rcpp::_) = Rcpp::CharacterVector::create(NA_STRING, NA_STRING);
    }
  }
  Rcpp::CharacterMatrix z_ptt = Rcpp::cbind(z, ppartition);

  double SS_ij = 0;
  Rcpp::List RnH;
  Rcpp::CharacterVector tmp;

  for (int i=0; i < M; i++) {
    m = ptts(i);
    arma::uvec idx_tmp = Rcpp::as<arma::uvec>(idx[m]);
    int n_ij = idx_tmp.n_elem;
    arma::uvec onevec = arma::ones<arma::uvec>(n_ij);
    arma::uvec idx_ij = idx_tmp - onevec;
    h = z_ptt(i,0);

    // prototype
    if (isin(m, ptts_proto)) {
      Rcpp::List RnH_list_byptt = RnH_list[m];
      RnH = RnH_list_byptt[h];
    } else {
      tmp = str_split(m, ",");
      std::string m1 = concatenate(tmp, "2", nd);
      Rcpp::List RnH_list_byptt = RnH_list[m1];
      RnH = RnH_list_byptt[h];
    }

    // children partitions
    Rcpp::CharacterVector cij = findchild(m, z_ptt, ptts);
    Rcpp::IntegerVector wherecij = match(cij, ptts);

    std::string mc;
    std::string hc;
    Rcpp::List cRnH;
    Rcpp::CharacterVector ctmp;

    // generate w
    arma::mat HRH = arma::zeros<arma::mat>(n_ij, n_ij);
    arma::vec HRw = arma::zeros<arma::vec>(n_ij);
    if (cij.length() > 0) {
      for (int ic=0; ic < cij.length(); ic++) {
        mc = cij(ic);
        arma::uvec idx_ctmp = Rcpp::as<arma::uvec>(idx[mc]);
        int n_cij = idx_ctmp.n_elem;
        arma::uvec onevec2 = arma::ones<arma::uvec>(n_cij);
        arma::uvec idx_cij = idx_ctmp - onevec2;
        hc = z_ptt(wherecij(ic)-1, 0);

        // prototype
        if (isin(mc, ptts_proto)) {
          Rcpp::List cRnH_list_byptt = RnH_list[mc];
          cRnH = cRnH_list_byptt[hc];
        } else {
          ctmp = str_split(mc, ",");
          std::string mc1 = concatenate(ctmp, "2", nd);
          Rcpp::List cRnH_list_byptt = RnH_list[mc1];
          cRnH = cRnH_list_byptt[hc];
        }

        // since the child partition can have other parent than m
        Rcpp::List cpptts_list_byptt = pptts_list[mc];
        Rcpp::List cpptts_list_byptt_bydir = cpptts_list_byptt[hc];
        arma::uvec idx_cptmp = Rcpp::as<arma::uvec>(cpptts_list_byptt_bydir["pidx"]);

        Rcpp::IntegerVector idx_mtmp = match(Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(idx_tmp)),
                                             Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(idx_cptmp)));
        int n_m = idx_mtmp.length();
        arma::uvec onevec3 = arma::ones<arma::uvec>(n_m);
        arma::uvec idx_m = Rcpp::as<arma::uvec>(idx_mtmp) - onevec3;

        arma::mat cRnH_H = cRnH["H"];
        arma::mat Hc = cRnH_H.cols(idx_m);
        arma::mat ctB = cRnH["tB"];
        arma::mat HR = Hc.t()*ctB.t()*ctB;
        HRH += HR*Hc;

        arma::uvec idx_nm = negidxing(idx_m,
                                      arma::linspace<arma::uvec>(0,idx_cptmp.n_elem-1,idx_cptmp.n_elem));

        arma::uvec idx_ncptmp = negidxing(idx_tmp, idx_cptmp);
        int n_ncp = idx_ncptmp.n_elem;
        arma::uvec onevec4 = arma::ones<arma::uvec>(n_ncp);
        arma::uvec idx_ncp = idx_ncptmp - onevec4;

        HRw += HR*(w(idx_cij) - cRnH_H.cols(idx_nm)*w(idx_ncp));
      }
    }

    // parent partitions
    Rcpp::List pptts_list_byptt = pptts_list[m];
    Rcpp::Nullable<Rcpp::List> pptts_list_byptt_bydir_ = pptts_list_byptt[h];
    arma::vec mu_ij;
    if (pptts_list_byptt_bydir_.isNotNull()) {
      // at least one parent
      Rcpp::List pptts_list_byptt_bydir(pptts_list_byptt_bydir_);
      arma::uvec idx_tmpp = Rcpp::as<arma::uvec>(pptts_list_byptt_bydir["pidx"]);
      int n_pij = idx_tmpp.n_elem;
      arma::uvec onevec5 = arma::ones<arma::uvec>(n_pij);
      arma::uvec idx_pij = idx_tmpp - onevec5;

      arma::mat RnH_H = RnH["H"];
      mu_ij = RnH_H*w(idx_pij);
    } else {
      // no parents
      mu_ij = arma::zeros<arma::vec>(n_ij);
    }

    arma::mat RnH_tB = RnH["tB"];
    arma::mat invR_ij = RnH_tB.t()*RnH_tB;
    HRH = symmatu(HRH);                                                         // had to force symmetry.

    arma::mat Vstar = arma::inv_sympd(invR_ij/sig_sq + arma::eye(n_ij, n_ij)/tau_sq + HRH/sig_sq);
    arma::vec mstar = (invR_ij*mu_ij)/sig_sq + res(idx_ij)/tau_sq + HRw/sig_sq;
    arma::vec w_ij = Vstar*mstar + arma::chol(Vstar, "lower")*arma::randn<arma::vec>(n_ij);
    w.elem(idx_ij) = w_ij;

    arma::vec res_ij = w_ij - mu_ij;
    arma::vec rtRr = res_ij.t()*invR_ij*res_ij;
    SS_ij += rtRr(0);
  }

  int n = w.n_elem;
  sig_sq = 1/arma::randg<double>(arma::distr_param(as + n/2, 1/(bs + SS_ij/2))); // scale parameter

  return Rcpp::List::create(
    Rcpp::Named("w") = w,
    Rcpp::Named("sig_sq") = sig_sq
  );
}

//////////////
// z_update //
//////////////
// [[Rcpp::export]]
Rcpp::CharacterVector zupdate(const arma::mat& pi_i, const arma::vec& w, const double& sig_sq,
                              const Rcpp::List& RnH_list, const Rcpp::CharacterVector& ptts,
                              const Rcpp::List& idx, const Rcpp::List& pptts_list,
                              const Rcpp::CharacterVector& ptts_proto, const double& nd,
                              const Rcpp::CharacterVector& directions) {

  std::string m;
  std::string h;
  int M = ptts.length();
  arma::mat ll(M, directions.length());
  Rcpp::CharacterVector out(M);

  for (int i=0; i < M; i++) {
    m = ptts(i);
    arma::uvec idx_tmp = Rcpp::as<arma::uvec>(idx[m]);
    int n_ij = idx_tmp.n_elem;
    arma::uvec onevec = arma::ones<arma::uvec>(n_ij);
    arma::uvec idx_ij = idx_tmp - onevec;
    Rcpp::List pptts_list_byptt = pptts_list[m];
    Rcpp::CharacterVector tmp;
    Rcpp::List RnH_wind;

    // prototype
    if (isin(m, ptts_proto)) {
      RnH_wind = RnH_list[m];
    } else {
      tmp = str_split(m, ",");
      std::string m1 = concatenate(tmp, "2", nd);
      RnH_wind = RnH_list[m1];
    }

    for (int j=0; j < directions.length(); j++) {
      h = directions(j);
      Rcpp::List RnH = RnH_wind[h];
      arma::mat tB = RnH["tB"];
      arma::mat invR = tB.t()*tB;
      arma::vec mu;

      Rcpp::Nullable<Rcpp::List> pptts_list_byptt_bydir_ = pptts_list_byptt[h];
      if (pptts_list_byptt_bydir_.isNotNull()) {
        Rcpp::List pptts_list_byptt_bydir(pptts_list_byptt_bydir_);
        arma::uvec idx_tmpp = Rcpp::as<arma::uvec>(pptts_list_byptt_bydir["pidx"]);
        int n_pij = idx_tmpp.n_elem;
        arma::uvec onevec2 = arma::ones<arma::uvec>(n_pij);
        arma::uvec idx_pij = idx_tmpp - onevec2;

        arma::mat H = RnH["H"];
        mu = w(idx_ij) - H*w(idx_pij);
      } else { mu = w(idx_ij); }

      arma::vec mtRm = mu.t()*invR*mu;
      ll(i,j) = log(pi_i(i,j)) + sum(log(diagvec(tB))) - 0.5*mtRm(0)/sig_sq;
    }

    arma::rowvec lprob = exp(ll.row(i) - max(ll.row(i)));
    Rcpp::NumericVector prob = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(lprob));
    Rcpp::CharacterVector samp = sample(directions, 1, false, prob);
    out(i) = samp(0);
  }
  return out;
}

////////////////
// psi_update //
////////////////
// [[Rcpp::export]]
Rcpp::List psiupdate(const arma::vec& w, const Rcpp::CharacterVector& z,
                     const double& sig_sq, arma::vec& theta, arma::mat& Sn,
                     const arma::mat& invS,
                     Rcpp::List& RnH_list,const arma::mat& coords,
                     const Rcpp::List& idx, const Rcpp::List& pptts_list,
                     const Rcpp::CharacterVector& ptts_proto,
                     const Rcpp::CharacterVector& ptts, const double& nd,
                     const double& iter, bool adaptive,
                     bool unifprior,
                     const Rcpp::CharacterVector& directions) {

  // propose freely from normal
  arma::vec U = arma::randn<arma::vec>(3);
  arma::vec theta_star = theta + Sn*U;

  double la, ua, lc, uc;
  arma::vec psi;
  arma::vec psi_star;
  if (unifprior) {
    la = invS(0,0);
    ua = invS(0,1);
    lc = invS(1,0);
    uc = invS(1,1);
    psi = {la + (ua-la)/(1+exp(-theta(0))),
           lc + (uc-lc)/(1+exp(-theta(1))),
           1/(1+exp(-theta(2)))};
    psi_star = {la + (ua-la)/(1+exp(-theta_star(0))),
                lc + (uc-lc)/(1+exp(-theta_star(1))),
                1/(1+exp(-theta_star(2)))};
  } else {
    psi = {exp(theta(0)), exp(theta(1)), 1/(1+exp(-theta(2)))};
    psi_star = {exp(theta_star(0)), exp(theta_star(1)), 1/(1+exp(-theta_star(2)))};
  }

  Rcpp::List RnH_star_list = createRnH(coords, idx, pptts_list, ptts_proto,
                                       psi_star(0), psi_star(1), psi_star(2),
                                       false, directions);

  double logr_star_ll = 0;
  double logr_ll = 0;
  std::string m;
  std::string h;
  for (int i=0; i < ptts.length(); i++) {
    m = ptts(i);
    arma::uvec idx_tmp = Rcpp::as<arma::uvec>(idx[m]);
    int n_ij = idx_tmp.n_elem;
    arma::uvec onevec = arma::ones<arma::uvec>(n_ij);
    arma::uvec idx_ij = idx_tmp - onevec;
    h = z(i);

    Rcpp::List RnH;
    Rcpp::List RnH_star;
    Rcpp::CharacterVector tmp;
    // prototype
    if (isin(m, ptts_proto)) {
      Rcpp::List RnH_list_byptt = RnH_list[m];
      RnH = RnH_list_byptt[h];
      Rcpp::List RnH_star_list_byptt = RnH_star_list[m];
      RnH_star = RnH_star_list_byptt[h];
    } else {
      tmp = str_split(m, ",");
      std::string m1 = concatenate(tmp, "2", nd);
      Rcpp::List RnH_list_byptt = RnH_list[m1];
      RnH = RnH_list_byptt[h];
      Rcpp::List RnH_star_list_byptt = RnH_star_list[m1];
      RnH_star = RnH_star_list_byptt[h];
    }

    // parent partitions
    Rcpp::List pptts_list_byptt = pptts_list[m];
    Rcpp::Nullable<Rcpp::List> pptts_list_byptt_bydir_ = pptts_list_byptt[h];
    arma::vec mu_ij;
    arma::vec mu_star_ij;
    if (pptts_list_byptt_bydir_.isNotNull()) {
      // at least one parent
      Rcpp::List pptts_list_byptt_bydir(pptts_list_byptt_bydir_);
      arma::uvec idx_tmpp = Rcpp::as<arma::uvec>(pptts_list_byptt_bydir["pidx"]);
      int n_pij = idx_tmpp.n_elem;
      arma::uvec onevec2 = arma::ones<arma::uvec>(n_pij);
      arma::uvec idx_pij = idx_tmpp - onevec2;

      arma::mat RnH_H = RnH["H"];
      mu_ij = RnH_H*w(idx_pij);

      arma::mat RnH_star_H = RnH_star["H"];
      mu_star_ij = RnH_star_H*w(idx_pij);
    } else {
      // no parents
      mu_ij = arma::zeros<arma::vec>(n_ij);
      mu_star_ij = arma::zeros<arma::vec>(n_ij);
    }
    arma::vec res_ij = w(idx_ij) - mu_ij;
    arma::vec res_star_ij = w(idx_ij) - mu_star_ij;

    arma::mat RnH_tB = RnH["tB"];
    arma::mat invR_ij = RnH_tB.t()*RnH_tB;

    arma::mat RnH_star_tB = RnH_star["tB"];
    arma::mat invR_star_ij = RnH_star_tB.t()*RnH_star_tB;

    arma::vec rtRr = res_ij.t()*invR_ij*res_ij;
    arma::vec rtRr_star = res_star_ij.t()*invR_star_ij*res_star_ij;
    logr_ll += sum(log(diagvec(RnH_tB))) - 0.5*rtRr(0)/sig_sq;
    logr_star_ll += sum(log(diagvec(RnH_star_tB))) - 0.5*rtRr_star(0)/sig_sq;
  }

  // accept or reject
  double logr;
  if (unifprior) {
    // Jacobian
    logr = logr_star_ll +
      log(psi_star(0)-la) + log(ua-psi_star(0)) +
      log(psi_star(1)-lc) + log(uc-psi_star(1)) +
      log(psi_star(2)) + log(1-psi_star(2)) -
      (logr_ll +
      log(psi(0)-la) + log(ua-psi(0)) +
      log(psi(1)-lc) + log(uc-psi(1)) +
      log(psi(2)) + log(1-psi(2)));
  } else {
    // normal prior
    arma::vec thtSth = theta.t()*invS*theta;
    arma::vec thtSth_star = theta_star.t()*invS*theta_star;
    logr = -0.5*thtSth_star(0) + logr_star_ll + 0.5*thtSth(0) - logr_ll;
  }

  if (log(arma::randu<double>()) < logr) {
    psi = psi_star;
    theta = theta_star;
    RnH_list = RnH_star_list;
  }

  // adaptive - update S
  if (adaptive) {
    double etan = std::min(1.0, 3*pow(iter, -2/3));
    double an = std::min(1.0, exp(logr));
    double U_l2norm = sum(square(U));
    Sn = chol(Sn*(arma::eye(3, 3) + U*U.t()*etan*(an - 0.234)/U_l2norm)*Sn.t(), "lower");
  }

  return Rcpp::List::create(
    Rcpp::Named("psi") = psi,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("Sn") = Sn,
    Rcpp::Named("RnH_list") = RnH_list
  );
}

