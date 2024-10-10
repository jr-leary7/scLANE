#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                     const Eigen::Map<Eigen::MatrixXd> B,
                     int n_cores = 1) {
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C;
  C.noalias() = A * B;
  return Rcpp::wrap(C);
}
