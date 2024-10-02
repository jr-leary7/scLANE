#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                     Eigen::Map<Eigen::MatrixXd> B,
                     int n_cores){

  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}
