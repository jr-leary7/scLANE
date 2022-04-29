#' Calculate part of the score statistic for a GLM.
#'
#' @description A function that calculates parts of the score statistic for GLMs only (it is used for the full path for forward selection).
#' @name stat_out_score_glm_null
#' @param Y : the response variable.
#' @param B_null : model matrix under the null model.
#' @return \code{stat_out_score_glm_null} returns a list of values (mainly products of matrices) that make up the final score statistic calculation (required for another function).
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats fitted.values
#' @importFrom gamlss gamlss
#' @seealso \code{\link{stat_out}} and \code{\link{stat_out_score_gee_null}}

stat_out_score_glm_null <- function(Y = NULL, B_null = NULL) {
  # check inputs
  if (is.null(Y) | is.null(B_null)) stop("Some inputs to stat_out_score_glm_null() are missing.")
  # run function
  ests <- gamlss::gamlss(Y ~ B_null - 1, family = "NBI", trace = FALSE)
  mu.est <- as.matrix(stats::fitted.values(ests))
  V.est <- mu.est * (1 + mu.est * (exp(ests$sigma.coef)))
  VS.est_list <- (c(Y) - c(mu.est)) / V.est
  # A_list <- chol2inv(chol((t(B_null) %*% diag(c(mu.est^2 / V.est)) %*% B_null)))
  temp_prod <- eigenMapMatMult(A = t(B_null),
                               B = diag(c(mu.est^2 / V.est)),
                               n_cores = 1)
  temp_prod <- eigenMapMatMult(A = temp_prod,
                               B = B_null,
                               n_cores = 1)
  A_list <- chol2inv(chol(temp_prod))
  # B1_list <- t(B_null) %*% diag(c(mu.est^2 / V.est))
  B1_list <- eigenMapMatMult(A = t(B_null),
                             B = diag(c(mu.est^2 / V.est)),
                             n_cores = 1)
  res <- list(VS.est_list = VS.est_list,
              A_list = A_list,
              B1_list = B1_list,
              mu.est = mu.est,
              V.est = V.est)
  return(res)
}
