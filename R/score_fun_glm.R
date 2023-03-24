#' Given estimates from the null model fit and the design matrix for alternative model, find the score statistic (this is used for GLMs only).
#'
#' @name score_fun_glm
#' @param Y : the response variable.
#' @param VS.est_list : a product of matrices.
#' @param A_list : a product of matrices.
#' @param B1_list : a product of matrices.
#' @param mu.est : estimates of the fitted mean under the null model.
#' @param V.est : estimates of the fitted variance under the null model.
#' @param B1 : model matrix under the null model.
#' @param XA : model matrix under the alternative model.
#' @return \code{score_fun_glm} returns a calculated score statistic for the null and alternative model when fitting a GLM.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats lm.fit

score_fun_glm <- function(Y = NULL,
                          VS.est_list = NULL,
                          A_list = NULL,
                          B1_list = NULL,
                          mu.est = NULL,
                          V.est = NULL,
                          B1 = NULL,
                          XA = NULL) {
  # check inputs
  if (is.null(Y) | is.null(VS.est_list) | is.null(A_list) | is.null(B1_list) | is.null(mu.est) | is.null(V.est) | is.null(B1) | is.null(XA)) { stop("Some inputs to score_fun_glm() are missing.") }
  # generate score statistic
  reg <- try(stats::lm.fit(B1, Y), silent = TRUE)  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.
  if (inherits(reg, "try-error")) {
    score <- NA
  } else if (any(is.na(reg$coef))) {
    score <- NA
  } else {
    VS.est_i <- unlist(VS.est_list)
    A_list_i <- unlist(A_list)
    B1_list_i <- unlist(B1_list)
    # B_list_i <- B1_list_i %*% XA
    B_list_i <- eigenMapMatMult(A = B1_list_i,
                                B = XA,
                                n_cores = 1)
    # D_list_i <- t(XA) %*% (XA * c((mu.est^2 / V.est)))
    D_list_i <- eigenMapMatMult(A = t(XA),
                                B = (XA * c((mu.est^2 / V.est))),
                                n_cores = 1)
    # inv.XVX_22 <- (D_list_i - t(B_list_i) %*% A_list_i %*% B_list_i)
    temp_prod <- eigenMapMatMult(A = t(B_list_i),
                                 B = A_list_i,
                                 n_cores = 1)
    temp_prod <- eigenMapMatMult(A = temp_prod,
                                 B = B_list_i,
                                 n_cores = 1)
    inv.XVX_22 <- (D_list_i - temp_prod); rm(temp_prod)
    # B.est <- t(((mu.est)) * VS.est_i) %*% XA
    B.est <- eigenMapMatMult(A = t(((mu.est)) * VS.est_i),
                             B = XA,
                             n_cores = 1)
    # score <- (B.est) %*% chol2inv(chol(inv.XVX_22)) %*% t(B.est)
    temp_prod <- eigenMapMatMult(A = B.est,
                                 B = chol2inv(chol(inv.XVX_22)),
                                 n_cores = 1)
    score <- eigenMapMatMult(A = temp_prod,
                             B = t(B.est),
                             n_cores = 1); rm(temp_prod)
  }
  res <- list(score = score)
  return(res)
}
