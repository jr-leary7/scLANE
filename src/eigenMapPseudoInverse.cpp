#include <RcppEigen.h>

#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT
#define EIGEN_DONT_VECTORIZE
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_UNROLLING_LIMIT 10
#define EIGEN_NO_MALLOC

using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd eigenMapPseudoInverse(const Eigen::Map<Eigen::MatrixXd> A,
                                      double tolerance = 1e-6,
                                      int n_cores = 1) {
  if (n_cores > 1) {
    Eigen::setNbThreads(n_cores);
  }

  JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
  MatrixXd singularValues = svd.singularValues();
  MatrixXd singularValuesInv(A.cols(), A.rows());
  singularValuesInv.setZero();
  for (int i = 0; i < singularValues.size(); ++i) {
    if (singularValues(i) > tolerance) {
      singularValuesInv(i, i) = 1.0 / singularValues(i);
    }
  }
  return svd.matrixV() * singularValuesInv * svd.matrixU().transpose();
}
