#' Given estimates from the null model fit and the design matrix for alternative model, find the score statistic (this is used for GEEs only).
#'
#' @name score_fun_gee
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @author Jack Leary
#' @importFrom stats lm.fit
#' @importFrom Matrix chol chol2inv
#' @description Calculate the score statistic for a GEE model.
#' @param Y The response variable. Defaults to NULL.
#' @param N The number of clusters. Defaults to NULL.
#' @param n_vec A vector consisting of the cluster sizes for each cluster. Defaults to NULL.
#' @param VS.est_list A product of matrices. Defaults to NULL.
#' @param AWA.est_list A product of matrices. Defaults to NULL.
#' @param J2_list A product of matrices. Defaults to NULL.
#' @param Sigma2_list A product of matrices. Defaults to NULL.
#' @param J11.inv A product of matrices. Defaults to NULL.
#' @param JSigma11 A product of matrices. Defaults to NULL.
#' @param mu.est Estimates of the fitted mean under the null model. Defaults to NULL.
#' @param V.est Estimates of the fitted variance under the null model. Defaults to NULL.
#' @param B1 Design matrix under the null model. Defaults to NULL.
#' @param XA Design matrix under the alternative model. Defaults to NULL.
#' @return A calculated score statistic for the null and alternative model when fitting a GEE.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @seealso \code{\link{score_fun_glm}}

score_fun_gee <- function(Y = NULL,
                          N = NULL,
                          n_vec = NULL,
                          VS.est_list = NULL,
                          AWA.est_list = NULL,
                          J2_list = NULL,
                          Sigma2_list = NULL,
                          J11.inv = NULL,
                          JSigma11 = NULL,
                          mu.est = NULL,
                          V.est = NULL,
                          B1 = NULL,
                          XA = NULL) {
  # check inputs
  if (is.null(Y) || is.null(N) || is.null(n_vec) || is.null(VS.est_list) || is.null(AWA.est_list) || is.null(J2_list) || is.null(Sigma2_list) || is.null(J11.inv) || is.null(JSigma11) || is.null(mu.est) || is.null(V.est) || is.null(B1) || is.null(XA)) { stop("Some inputs to score_fun_gee() are missing.") }
  # generate score statistic
  reg <- try(stats::lm.fit(B1, Y), silent = TRUE)  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.
  if (inherits(reg, "try-error")) {
    score <- NA
  } else if (any(is.na(reg$coef))) {
    score <- NA
  } else {
    p <- ncol(XA)
    p1 <- ncol(B1) - ncol(XA)
    n_vec1 <- c(0, n_vec)
    B.est <- matrix(0, nrow = p, ncol = 1)
    Sigma22 <- matrix(0, nrow = p, ncol = p)
    J21 <- Sigma21 <- matrix(0, nrow = p, ncol = p1)

    for (i in seq(N)) {
      k <- sum(n_vec[seq(i)])
      VS.est_i <- VS.est_list[[i]]
      AWA.est_i <- AWA.est_list[[i]]
      J2_i <- J2_list[[i]]
      Sigma2_i <- Sigma2_list[[i]]
      D.est_i <- eigenMapMatMult(A = diag((mu.est[(sum(n_vec1[seq(i)]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i]),
                                 B = XA[(sum(n_vec1[seq(i)]) + 1):k, ],
                                 n_cores = 1)
      D_est_i_transpose <- t(D.est_i)
      J21 <- J21 + eigenMapMatMult(A = D_est_i_transpose,
                                   B = t(J2_i),
                                   n_cores = 1)
      Sigma21 <- Sigma21 + eigenMapMatMult(A = D_est_i_transpose,
                                           B = t(Sigma2_i),
                                           n_cores = 1)
      B.est <- B.est + eigenMapMatMult(A = D_est_i_transpose,
                                       B = VS.est_i,
                                       n_cores = 1)
      temp_prod <- eigenMapMatMult(A = D_est_i_transpose,
                                   B = AWA.est_i,
                                   n_cores = 1)
      Sigma22 <- Sigma22 + eigenMapMatMult(A = temp_prod,
                                           B = D.est_i,
                                           n_cores = 1)
    }
    temp_prod_1 <- eigenMapMatMult(A = J21,
                                   B = J11.inv,
                                   n_cores = 1)
    temp_prod_1 <- eigenMapMatMult(A = temp_prod_1,
                                   B = t(Sigma21),
                                   n_cores = 1)
    temp_prod_2 <- eigenMapMatMult(A = Sigma21,
                                   B = J11.inv,
                                   n_cores = 1)
    J21_transpose <- t(J21)
    temp_prod_2 <- eigenMapMatMult(A = temp_prod_2,
                                   B = J21_transpose,
                                   n_cores = 1)
    temp_prod_3 <- eigenMapMatMult(A = J21,
                                   B = JSigma11,
                                   n_cores = 1)
    temp_prod_3 <- eigenMapMatMult(A = temp_prod_3,
                                   B = J21_transpose,
                                   n_cores = 1)
    Sigma <- Sigma22 - temp_prod_1 - temp_prod_2 + temp_prod_3
    temp_prod <- eigenMapMatMult(A = t(B.est),
                                 B = Matrix::chol2inv(Matrix::chol(Sigma)),
                                 n_cores = 1)
    score <- eigenMapMatMult(A = temp_prod,
                             B = B.est,
                             n_cores = 1)
  }
  res <- list(score = score)
  return(res)
}
