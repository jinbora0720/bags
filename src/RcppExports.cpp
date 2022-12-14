// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bag_predict
Rcpp::List bag_predict(const arma::mat& coords_tr, const arma::mat& coords_pred, const Rcpp::List& ptts_pred, const Rcpp::List& idx_pred, const Rcpp::List& pptts_pred_list, const arma::mat& w_save, const arma::mat& psi_save, const arma::vec& sigsq_save, const arma::vec& tausq_save, int num_threads, bool verbose);
RcppExport SEXP _bags_bag_predict(SEXP coords_trSEXP, SEXP coords_predSEXP, SEXP ptts_predSEXP, SEXP idx_predSEXP, SEXP pptts_pred_listSEXP, SEXP w_saveSEXP, SEXP psi_saveSEXP, SEXP sigsq_saveSEXP, SEXP tausq_saveSEXP, SEXP num_threadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_tr(coords_trSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords_pred(coords_predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type ptts_pred(ptts_predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx_pred(idx_predSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pptts_pred_list(pptts_pred_listSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_save(w_saveSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi_save(psi_saveSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigsq_save(sigsq_saveSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tausq_save(tausq_saveSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(bag_predict(coords_tr, coords_pred, ptts_pred, idx_pred, pptts_pred_list, w_save, psi_save, sigsq_save, tausq_save, num_threads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// createRnH
Rcpp::List createRnH(const arma::mat& coords, const Rcpp::List& idx, const Rcpp::List& pptts_list, const Rcpp::CharacterVector& ptts_proto, const double& a, const double& c, const double& kappa, bool produce_R, const Rcpp::CharacterVector& directions);
RcppExport SEXP _bags_createRnH(SEXP coordsSEXP, SEXP idxSEXP, SEXP pptts_listSEXP, SEXP ptts_protoSEXP, SEXP aSEXP, SEXP cSEXP, SEXP kappaSEXP, SEXP produce_RSEXP, SEXP directionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pptts_list(pptts_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts_proto(ptts_protoSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< bool >::type produce_R(produce_RSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type directions(directionsSEXP);
    rcpp_result_gen = Rcpp::wrap(createRnH(coords, idx, pptts_list, ptts_proto, a, c, kappa, produce_R, directions));
    return rcpp_result_gen;
END_RCPP
}
// wsigupdate
Rcpp::List wsigupdate(const arma::vec& res, arma::vec& w, double& sig_sq, const double& tau_sq, const Rcpp::CharacterVector& z, const Rcpp::List& RnH_list, const Rcpp::CharacterVector& ptts, const Rcpp::List& idx, const Rcpp::List& pptts_list, const Rcpp::CharacterVector& ptts_proto, const double& nd, const double& as, const double& bs);
RcppExport SEXP _bags_wsigupdate(SEXP resSEXP, SEXP wSEXP, SEXP sig_sqSEXP, SEXP tau_sqSEXP, SEXP zSEXP, SEXP RnH_listSEXP, SEXP pttsSEXP, SEXP idxSEXP, SEXP pptts_listSEXP, SEXP ptts_protoSEXP, SEXP ndSEXP, SEXP asSEXP, SEXP bsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double& >::type sig_sq(sig_sqSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_sq(tau_sqSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type RnH_list(RnH_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts(pttsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pptts_list(pptts_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts_proto(ptts_protoSEXP);
    Rcpp::traits::input_parameter< const double& >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< const double& >::type as(asSEXP);
    Rcpp::traits::input_parameter< const double& >::type bs(bsSEXP);
    rcpp_result_gen = Rcpp::wrap(wsigupdate(res, w, sig_sq, tau_sq, z, RnH_list, ptts, idx, pptts_list, ptts_proto, nd, as, bs));
    return rcpp_result_gen;
END_RCPP
}
// zupdate
Rcpp::CharacterVector zupdate(const arma::mat& pi_i, const arma::vec& w, const double& sig_sq, const Rcpp::List& RnH_list, const Rcpp::CharacterVector& ptts, const Rcpp::List& idx, const Rcpp::List& pptts_list, const Rcpp::CharacterVector& ptts_proto, const double& nd, const Rcpp::CharacterVector& directions);
RcppExport SEXP _bags_zupdate(SEXP pi_iSEXP, SEXP wSEXP, SEXP sig_sqSEXP, SEXP RnH_listSEXP, SEXP pttsSEXP, SEXP idxSEXP, SEXP pptts_listSEXP, SEXP ptts_protoSEXP, SEXP ndSEXP, SEXP directionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type pi_i(pi_iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_sq(sig_sqSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type RnH_list(RnH_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts(pttsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pptts_list(pptts_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts_proto(ptts_protoSEXP);
    Rcpp::traits::input_parameter< const double& >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type directions(directionsSEXP);
    rcpp_result_gen = Rcpp::wrap(zupdate(pi_i, w, sig_sq, RnH_list, ptts, idx, pptts_list, ptts_proto, nd, directions));
    return rcpp_result_gen;
END_RCPP
}
// psiupdate
Rcpp::List psiupdate(const arma::vec& w, const Rcpp::CharacterVector& z, const double& sig_sq, arma::vec& theta, arma::mat& Sn, const arma::mat& invS, Rcpp::List& RnH_list, const arma::mat& coords, const Rcpp::List& idx, const Rcpp::List& pptts_list, const Rcpp::CharacterVector& ptts_proto, const Rcpp::CharacterVector& ptts, const double& nd, const double& iter, bool adaptive, bool unifprior, const Rcpp::CharacterVector& directions);
RcppExport SEXP _bags_psiupdate(SEXP wSEXP, SEXP zSEXP, SEXP sig_sqSEXP, SEXP thetaSEXP, SEXP SnSEXP, SEXP invSSEXP, SEXP RnH_listSEXP, SEXP coordsSEXP, SEXP idxSEXP, SEXP pptts_listSEXP, SEXP ptts_protoSEXP, SEXP pttsSEXP, SEXP ndSEXP, SEXP iterSEXP, SEXP adaptiveSEXP, SEXP unifpriorSEXP, SEXP directionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type sig_sq(sig_sqSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sn(SnSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type RnH_list(RnH_listSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pptts_list(pptts_listSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts_proto(ptts_protoSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts(pttsSEXP);
    Rcpp::traits::input_parameter< const double& >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< const double& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< bool >::type unifprior(unifpriorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type directions(directionsSEXP);
    rcpp_result_gen = Rcpp::wrap(psiupdate(w, z, sig_sq, theta, Sn, invS, RnH_list, coords, idx, pptts_list, ptts_proto, ptts, nd, iter, adaptive, unifprior, directions));
    return rcpp_result_gen;
END_RCPP
}
// isin
bool isin(const std::string& m, const Rcpp::CharacterVector& ptts_proto);
RcppExport SEXP _bags_isin(SEXP mSEXP, SEXP ptts_protoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts_proto(ptts_protoSEXP);
    rcpp_result_gen = Rcpp::wrap(isin(m, ptts_proto));
    return rcpp_result_gen;
END_RCPP
}
// str_split
std::vector<std::string> str_split(const std::string& str, const char* delim);
RcppExport SEXP _bags_str_split(SEXP strSEXP, SEXP delimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type str(strSEXP);
    Rcpp::traits::input_parameter< const char* >::type delim(delimSEXP);
    rcpp_result_gen = Rcpp::wrap(str_split(str, delim));
    return rcpp_result_gen;
END_RCPP
}
// concatenate
std::string concatenate(const Rcpp::CharacterVector& x, const std::string& z, const double& nd);
RcppExport SEXP _bags_concatenate(SEXP xSEXP, SEXP zSEXP, SEXP ndSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type nd(ndSEXP);
    rcpp_result_gen = Rcpp::wrap(concatenate(x, z, nd));
    return rcpp_result_gen;
END_RCPP
}
// findchild
Rcpp::CharacterVector findchild(const std::string& m, const Rcpp::CharacterMatrix& z_ptt, const Rcpp::CharacterVector& ptts);
RcppExport SEXP _bags_findchild(SEXP mSEXP, SEXP z_pttSEXP, SEXP pttsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterMatrix& >::type z_ptt(z_pttSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type ptts(pttsSEXP);
    rcpp_result_gen = Rcpp::wrap(findchild(m, z_ptt, ptts));
    return rcpp_result_gen;
END_RCPP
}
// negidxing
arma::uvec negidxing(const arma::uvec& a, const arma::uvec& A);
RcppExport SEXP _bags_negidxing(SEXP aSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negidxing(a, A));
    return rcpp_result_gen;
END_RCPP
}
// distmat_Euclidean
arma::mat distmat_Euclidean(const arma::mat& A);
RcppExport SEXP _bags_distmat_Euclidean(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(distmat_Euclidean(A));
    return rcpp_result_gen;
END_RCPP
}
// distmat_abs
arma::mat distmat_abs(const arma::mat& A);
RcppExport SEXP _bags_distmat_abs(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(distmat_abs(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bags_bag_predict", (DL_FUNC) &_bags_bag_predict, 11},
    {"_bags_createRnH", (DL_FUNC) &_bags_createRnH, 9},
    {"_bags_wsigupdate", (DL_FUNC) &_bags_wsigupdate, 13},
    {"_bags_zupdate", (DL_FUNC) &_bags_zupdate, 10},
    {"_bags_psiupdate", (DL_FUNC) &_bags_psiupdate, 17},
    {"_bags_isin", (DL_FUNC) &_bags_isin, 2},
    {"_bags_str_split", (DL_FUNC) &_bags_str_split, 2},
    {"_bags_concatenate", (DL_FUNC) &_bags_concatenate, 3},
    {"_bags_findchild", (DL_FUNC) &_bags_findchild, 3},
    {"_bags_negidxing", (DL_FUNC) &_bags_negidxing, 2},
    {"_bags_distmat_Euclidean", (DL_FUNC) &_bags_distmat_Euclidean, 1},
    {"_bags_distmat_abs", (DL_FUNC) &_bags_distmat_abs, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_bags(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
