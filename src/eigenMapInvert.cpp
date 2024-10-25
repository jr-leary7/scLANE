#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd eigenMapMatrixInvert(const Eigen::Map<Eigen::MatrixXd> A, int n_cores = 1) {
  if (n_cores > 1) {
    Eigen::setNbThreads(n_cores);
  }
  Eigen::MatrixXd B = A.llt().solve(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
  return B;
}
