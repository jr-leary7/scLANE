#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

SEXP eigenMapMatrixInvert(const Eigen::Map<Eigen::MatrixXd> A, int n_cores){
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd B = A.llt().solve(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
  return Rcpp::wrap(B);
}
