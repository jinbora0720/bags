// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
bool isin(const std::string& m, const Rcpp::CharacterVector& ptts_proto) {

  int n = ptts_proto.length();
  bool out = false;
  for (int i=0; i < n; i++) {
    std::string m1 = Rcpp::as<std::string>(ptts_proto(i));
    if (m == m1) { out = true; }
  }

  return out;
}

// [[Rcpp::export]]
std::vector<std::string> str_split(const std::string& str, const char* delim) {
  std::vector<std::string> strings;
  size_t start;
  size_t end = 0;
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    end = str.find(delim, start);
    strings.push_back(str.substr(start, end - start));
  }
  return strings;
}

// [[Rcpp::export]]
std::string concatenate(const Rcpp::CharacterVector& x,
                        const std::string& z, const double& nd) {
  std::string x1;
  x1 = x(0);
  std::string x2;
  x2 = x(1);
  std::string z1;
  for (int i=0; i < nd-1; i++) {
    z1 += "0";
  }
  z1 += z;
  return x1 + "," + x2 + "," + z1;
}

// [[Rcpp::export]]
Rcpp::CharacterVector findchild(const std::string& m, const Rcpp::CharacterMatrix& z_ptt,
                                const Rcpp::CharacterVector& ptts) {

  int n = z_ptt.nrow();
  Rcpp::CharacterVector out;
  for (int i=0; i < n; i++) {
    std::string m1 = Rcpp::as<std::string>(z_ptt(i,1));
    std::string m2 = Rcpp::as<std::string>(z_ptt(i,2));
    if (m == m1 || m == m2) {
      std::string name = Rcpp::as<std::string>(ptts(i));
      out.push_back(name);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::uvec negidxing(const arma::uvec& a, const arma::uvec& A) {
  // remove a from A
  int N = A.n_elem;
  int n = a.n_elem;
  Rcpp::IntegerVector out;

  for (int i=0; i < N; i++) {
    int Ai = A(i);
    double test = 0;
    for (int j=0; j < n; j++) {
      int aj = a(j);
      if (Ai == aj) { test = 1; }
    }
    if (test == 0) { out.push_back(Ai); }
  }

  return Rcpp::as<arma::uvec>(out);
}
