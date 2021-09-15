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
#' @param ... : further arguments passed to or from other methods.
#' @return \code{score_fun_glm} returns a calculated score statistic for the null and alternative model when fitting a GLM.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @export
#' @seealso \code{\link{score_fun_gee}}

score_fun_glm <- function(Y = NULL,
                          VS.est_list = NULL,
                          A_list = NULL,
                          B1_list = NULL,
                          mu.est = NULL,
                          V.est = NULL,
                          B1 = NULL,
                          XA = NULL,
                          ...) {
  # check inputs
  if (any(sapply(c(Y, VS.est_list, A_list, B1_list, mu.est, V.est, B1, XA), is.null))) stop("Some inputs are missing from score_fun_glm().")
  # generate score statistic
  reg <- try(stats::lm.fit(B1, Y, ...), silent = TRUE)  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.
  if (class(reg)[1] == "try-error") {
    score <- NA
  } else if (any(is.na(reg$coef))) {
    score <- NA
  } else {
    N <- length(Y)
    VS.est_i <- unlist(VS.est_list)
    A_list_i <- unlist(A_list)
    B1_list_i <- unlist(B1_list)
    B_list_i <- B1_list_i %*% XA
    D_list_i <- t(XA) %*% (XA * c((mu.est^2 / V.est)))
    inv.XVX_22 <- (D_list_i - t(B_list_i) %*% A_list_i %*% B_list_i)
    B.est <- t(((mu.est)) * VS.est_i) %*% XA
    if (N < 1500) {
      score <- (B.est) %*% MASS::ginv(inv.XVX_22) %*% t(B.est)
    } else {
      score <- (B.est) %*% chol2inv(chol(inv.XVX_22)) %*% t(B.est)
    }
  }
  res <- list(score = score)
  return(res)
}
