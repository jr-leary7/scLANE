// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// eigenMapMatrixInvert
Eigen::MatrixXd eigenMapMatrixInvert(const Eigen::Map<Eigen::MatrixXd> A, int n_cores);
RcppExport SEXP _scLANE_eigenMapMatrixInvert(SEXP ASEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapMatrixInvert(A, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// eigenMapMatMult
Eigen::MatrixXd eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B, int n_cores);
RcppExport SEXP _scLANE_eigenMapMatMult(SEXP ASEXP, SEXP BSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapMatMult(A, B, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// eigenMapPseudoInverse
Eigen::MatrixXd eigenMapPseudoInverse(const Eigen::Map<Eigen::MatrixXd> A, double tolerance, int n_cores);
RcppExport SEXP _scLANE_eigenMapPseudoInverse(SEXP ASEXP, SEXP toleranceSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapPseudoInverse(A, tolerance, n_cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scLANE_eigenMapMatrixInvert", (DL_FUNC) &_scLANE_eigenMapMatrixInvert, 2},
    {"_scLANE_eigenMapMatMult", (DL_FUNC) &_scLANE_eigenMapMatMult, 3},
    {"_scLANE_eigenMapPseudoInverse", (DL_FUNC) &_scLANE_eigenMapPseudoInverse, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scLANE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
