#' Given estimates from the null model fit and the design matrix for alternative model, find the score statistic (this is used for GEEs only).
#'
#' @name score_fun_gee
#' @param Y : the response variable.
#' @param N : the number of clusters.
#' @param n_vec : a vector consisting of the cluster sizes for each cluster.
#' @param VS.est_list : a product of matrices.
#' @param AWA.est_list : a product of matrices.
#' @param J2_list : a product of matrices.
#' @param Sigma2_list : a product of matrices.
#' @param J11.inv : a product of matrices.
#' @param JSigma11 : a product of matrices.
#' @param mu.est : estimates of the fitted mean under the null model.
#' @param V.est : estimates of the fitted variance under the null model.
#' @param B1 : model matrix under the null model.
#' @param XA : model matrix under the alternative model.
#' @return \code{score_fun_gee} returns a calculated score statistic for the null and alternative model when fitting a GEE.
#' @author Jakub Stoklosa and David I. Warton
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats lm.fit
#' @importFrom MASS ginv
#' @export
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
  if (is.null(Y) | is.null(N) | is.null(n_vec) | is.null(VS.est_list) | is.null(AWA.est_list) | is.null(J2_list) | is.null(Sigma2_list) | is.null(J11.inv) | is.null(JSigma11) | is.null(mu.est) | is.null(V.est) | is.null(B1) | is.null(XA)) { stop("Some inputs to score_fun_gee() are missing.") }
  # generate score statistic
  reg <- try(stats::lm.fit(B1, Y), silent = TRUE)  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.
  if (class(reg)[1] == "try-error") {
    score <- NA
  } else if (any(is.na(reg$coef))) {
    score <- NA
  } else {
    p <- ncol(XA)
    p1 <- ncol(B1) - ncol(XA)
    n_vec1 <- c(0, n_vec)
    B.est <- matrix(0, nrow = p, ncol = 1)
    Sigma22 <- matrix(0, nrow = p, ncol = p)
    J21 <- matrix(0, nrow = p, ncol = p1)
    Sigma21 <- matrix(0, nrow = p, ncol = p1)

    for (i in seq(N)) {
      k <- sum(n_vec[1:i])
      VS.est_i <- VS.est_list[[i]]
      AWA.est_i <- AWA.est_list[[i]]
      J2_i <- J2_list[[i]]
      Sigma2_i <- Sigma2_list[[i]]
      D.est_i <- diag((mu.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i]) %*% XA[(sum(n_vec1[1:i]) + 1):k, ]
      J21 <- J21 + t(D.est_i) %*% t(J2_i)
      Sigma21 <- Sigma21 + t(D.est_i) %*% t(Sigma2_i)
      B.est <- B.est + t(D.est_i) %*% VS.est_i
      Sigma22 <- Sigma22 + t(D.est_i) %*% AWA.est_i %*% (D.est_i)
    }
    Sigma <- Sigma22 - (J21 %*% J11.inv) %*% t(Sigma21) - (Sigma21 %*% J11.inv) %*% t(J21) + J21 %*% JSigma11 %*% t(J21)
    score <- t(B.est) %*% MASS::ginv(Sigma) %*% B.est
  }
  res <- list(score = score)
  return(res)
}