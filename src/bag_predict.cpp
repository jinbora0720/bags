#ifdef _OPENMP
#include <omp.h>
#endif

#include "matrix_distance.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::export]]
Rcpp::List bag_predict(const arma::mat& coords_tr,
                       const arma::mat& coords_pred,
                       const Rcpp::List& ptts_pred,
                       const Rcpp::List& idx_pred,
                       const Rcpp::List& pptts_pred_list,
                       const arma::mat& w_save,
                       const arma::mat& psi_save,
                       const arma::vec& sigsq_save,
                       const arma::vec& tausq_save,
                       int num_threads = 4,
                       bool verbose = false) {

#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif

  int n = coords_pred.n_rows;
  int save = psi_save.n_cols;
  int L = ptts_pred.length();
  arma::uvec ColInd = {0, 1, 2};
  arma::mat w_pred(n, save);
  arma::mat y_pred(n, save);

  Progress p(L, verbose);

  for(int i=0; i < L; i++){
    if (verbose) {
      p.increment(); // update progress
    }

    std::string l = ptts_pred(i);
    arma::uvec idx_l = Rcpp::as<arma::uvec>(idx_pred[l]);
    int n_l = idx_l.n_elem;
    arma::uvec onevec = arma::ones<arma::uvec>(n_l);
    idx_l -= onevec;

    Rcpp::List pptts_pred_list_byptt = pptts_pred_list[l];
    arma::uvec idx_pl = Rcpp::as<arma::uvec>(pptts_pred_list_byptt["pidx"]);
    int n_pl = idx_pl.n_elem;
    arma::uvec onevec2 = arma::ones<arma::uvec>(n_pl);
    idx_pl -= onevec2;

    arma::mat coords_tmp = join_cols(coords_pred.submat(idx_l, ColInd),
                               coords_tr.submat(idx_pl, ColInd));
    arma::mat spatdist = distmat_Euclidean(coords_tmp.submat(0, 0, n_l + n_pl -1, 1));
    arma::mat timedist = distmat_abs(coords_tmp.submat(0, 2, n_l + n_pl -1, 2));

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int s=0; s < save; s++) {
      arma::vec psi = psi_save.col(s);
      arma::vec w = w_save.col(s);
      double sigsq = sigsq_save(s);
      double tausq = tausq_save(s);
      arma::mat aup1 = psi(0)*timedist+1;
      arma::mat Cor = exp(-psi(1)*spatdist/pow(aup1, psi(2)/2))/aup1;
      arma::mat Cipi = Cor.submat(0, n_l, n_l-1, n_l+n_pl-1);
      arma::mat CinvC = Cipi*
        inv_sympd(Cor.submat(n_l, n_l, n_l+n_pl-1, n_l+n_pl-1));
      arma::vec CinvCC_diag = diagvec(CinvC*Cipi.t());
      arma::vec w_mean = CinvC *w(idx_pl);

      for (int j=0; j < n_l; j++) {
        int row = idx_l(j);
        double sqrtR = pow(1 - CinvCC_diag(j), 0.5);
        double w_lj = w_mean(j) +
          pow(sigsq, .5)*sqrtR*arma::conv_to<double>::from(arma::randn(1));
        w_pred(row, s) = w_lj;
        y_pred(row, s) = w_lj + pow(tausq, .5)*arma::conv_to<double>::from(arma::randn(1));
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("w_pred_save") = w_pred,
    Rcpp::Named("y_pred_save") = y_pred
  );
}
