
#include "matrix_distance.h"

// [[Rcpp::export]]
arma::mat distmat_Euclidean(const arma::mat& A) {

  int n = A.n_rows;
  arma::mat out = arma::zeros<arma::mat>(n, n);

  for (int i=0; i < n; i++) {
    for (int j=i+1; j < n; j++) {
      out(i, j) = norm((A.row(i) - A.row(j)));
    }
  }
  out += out.t();

  return out;
}

// [[Rcpp::export]]
arma::mat distmat_abs(const arma::mat& A) {

  int n = A.n_rows;
  arma::mat out = arma::zeros<arma::mat>(n, n);

  for (int i=0; i < n; i++) {
    for (int j=i+1; j < n; j++) {
      double diff = A(i,0) - A(j,0);
      out(i, j) = std::abs(diff);
    }
  }
  out += out.t();

  return out;
}
