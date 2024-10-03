#' Given estimates from the null model fit and the design matrix for alternative model, find the score statistic (this is used for GLMs only).
#'
#' @name score_fun_glm
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @author Jack Leary
#' @importFrom stats lm.fit
#' @importFrom MASS ginv
#' @description Calculate the score statistic for a GLM model.
#' @param Y The response variable. Defaults to NULL.
#' @param VS.est_list A product of matrices. Defaults to NULL.
#' @param A_list A product of matrices. Defaults to NULL.
#' @param B1_list A product of matrices. Defaults to NULL.
#' @param mu.est Estimates of the fitted mean under the null model. Defaults to NULL.
#' @param V.est Estimates of the fitted variance under the null model. Defaults to NULL.
#' @param B1 Design matrix under the null model. Defaults to NULL.
#' @param XA Design matrix under the alternative model. Defaults to NULL.
#' @return A calculated score statistic for the null and alternative model when fitting a GLM.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.

score_fun_glm <- function(Y = NULL,
                          VS.est_list = NULL,
                          A_list = NULL,
                          B1_list = NULL,
                          mu.est = NULL,
                          V.est = NULL,
                          B1 = NULL,
                          XA = NULL) {
  # check inputs
  if (is.null(Y) || is.null(VS.est_list) || is.null(A_list) || is.null(B1_list) || is.null(mu.est) || is.null(V.est) || is.null(B1) || is.null(XA)) { stop("Some inputs to score_fun_glm() are missing.") }
  # generate score statistic
  reg <- try({ stats::lm.fit(B1, Y) }, silent = TRUE)  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.
  if (inherits(reg, "try-error") || any(is.na(reg$coefficients))) {
    score <- NA_real_
  } else {
    VS.est_i <- unlist(VS.est_list)
    A_list_i <- unlist(A_list)
    B1_list_i <- unlist(B1_list)
    B_list_i <- eigenMapMatMult(A = B1_list_i,
                                B = XA,
                                n_cores = 1)
    D_list_i <- eigenMapMatMult(A = t(XA),
                                B = (XA * c(mu.est^2 / V.est)),
                                n_cores = 1)
    temp_prod <- eigenMapMatMult(A = t(B_list_i),
                                 B = A_list_i,
                                 n_cores = 1)
    temp_prod <- eigenMapMatMult(A = temp_prod,
                                 B = B_list_i,
                                 n_cores = 1)
    inv.XVX_22 <- D_list_i - temp_prod
    B.est <- eigenMapMatMult(A = t(mu.est * VS.est_i),
                             B = XA,
                             n_cores = 1)
    XVX_22 <- try({ eigenMapMatrixInvert(inv.XVX_22, n_cores = 1) }, silent = TRUE)
    if (inherits(XVX_22, "try-error")) {
      XVX_22 <- MASS::ginv(inv.XVX_22)
    }
    temp_prod <- eigenMapMatMult(A = B.est,
                                 B = XVX_22,
                                 n_cores = 1)
    score <- eigenMapMatMult(A = temp_prod,
                             B = t(B.est),
                             n_cores = 1)
  }
  res <- list(score = score)
  return(res)
}
