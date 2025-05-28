#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd eigenMapPseudoInverse(const Eigen::Map<Eigen::MatrixXd> A,
                                      double tolerance = 1e-6,
                                      int n_cores = 1) {
  if (n_cores > 1) {
    Eigen::setNbThreads(n_cores);
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd singularValues = svd.singularValues();
  Eigen::MatrixXd singularValuesInv(A.cols(), A.rows());
  singularValuesInv.setZero();
  for (int i = 0; i < singularValues.size(); ++i) {
    if (singularValues(i) > tolerance) {
      singularValuesInv(i, i) = 1.0 / singularValues(i);
    }
  }
  return svd.matrixV() * singularValuesInv * svd.matrixU().transpose();
}
