#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A,
                                const Eigen::Map<Eigen::MatrixXd> B,
                                int n_cores = 1) {
  if (n_cores > 1) {
    Eigen::setNbThreads(n_cores);
  }
  Eigen::MatrixXd C;
  C.noalias() = A * B;
  return C;
}
