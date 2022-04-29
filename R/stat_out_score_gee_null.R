#' Calculate part of the score statistic for a GEE
#'
#' @description A function that calculates parts of the score statistic for GEEs only (it is used for the full path for forward selection).
#' @name stat_out_score_gee_null
#' @param Y : the response variable.
#' @param B_null : model matrix under the null model.
#' @param id.vec If \code{is.gee = TRUE}, must be a vector of ID values for the observations. Defaults to NULL.
#' @param cor.structure If \code{is.gee = TRUE}, must be a string specifying the desired correlation structure for the NB GEE. Defaults to NULL.
#' @return \code{stat_out_score_glm_null} returns a list of values (mainly products of matrices) that make up the final score statistic calculation (required for another function).
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom geeM geem
#' @importFrom gamlss gamlss random
#' @importFrom stats fitted.values
#' @seealso \code{\link{stat_out}} and \code{\link{stat_out_score_glm_null}}

stat_out_score_gee_null <- function(Y = NULL,
                                    B_null = NULL,
                                    id.vec = NULL,
                                    cor.structure = NULL) {
  # check inputs
  if (is.null(Y) | is.null(B_null) | is.null(id.vec) | is.null(cor.structure)) { stop("All arguments to stat_out_score_gee_null() must be non-NULL.") }
  cor.structure <- tolower(cor.structure)
  if (!cor.structure %in% c("independence", "exchangeable", "ar1")) { stop("cor.structure in stat_out_score_gee_null() must be a known type.") }
  n_vec <- as.numeric(table(id.vec))
  N <- length(unique(id.vec))
  # fit null NB GEE
  ests <- geeM::geem(Y ~ 1,
                     id = id.vec,
                     corstr = cor.structure,
                     family = MASS::negative.binomial(1),
                     sandwich = TRUE,
                     maxit = 1)
  ests_gam <- gamlss::gamlss(Y ~ gamlss::random(as.factor(id.vec)),
                             family = "NBI",
                             trace = FALSE)
  alpha_est <- ests$alpha
  mu_est <- as.matrix(stats::fitted.values(ests))
  V_est <- mu_est * (1 + mu_est * (exp(ests_gam$sigma.coefficients)))  # NB var = mu + mu^2 * sigma

  p <- ncol(B_null)
  n_vec1 <- c(0, n_vec)
  VS_est_list <- list()
  AWA_est_list <- list()
  J2_list <- list()
  Sigma2_list <- list()
  J11 <- matrix(0, nrow = p, ncol = p)
  Sigma11 <- matrix(0, nrow = p, ncol = p)

  for (i in seq(N)) {
    k <- sum(n_vec[1:i])
    # set up working correlation matrix structure
    if (cor.structure == "independence") {
      R_alpha <- diag(1, nrow = n_vec[i], ncol = n_vec[i])
    } else if (cor.structure == "ar1") {
      R_alpha <- matrix(c(rep(alpha_est, n_vec[i]*n_vec[i])), ncol = n_vec[i]) + diag(c(1 - alpha_est), ncol = n_vec[i], nrow = n_vec[i])
    } else if (cor.structure == "exchangeable") {
      R_alpha <- alpha_est^outer(1:n_vec[i], 1:n_vec[i], function(x, y) abs(x - y))
    } else {
      stop("Currently unsupported correlation structure.")
    }

    V_est_i <- diag(sqrt(V_est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i]) %*% R_alpha %*% diag(sqrt(V_est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])
    V_est_i_inv <- chol2inv(chol(V_est_i))
    S_est_i <- c(t(Y))[(sum(n_vec1[1:i]) + 1):k] - mu_est[(sum(n_vec1[1:i]) + 1):k]
    AWA_est_i <- V_est_i_inv %*% (S_est_i %*% t(S_est_i)) %*% V_est_i_inv
    D_est_i <- diag((mu_est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i]) %*% B_null[(sum(n_vec1[1:i]) + 1):k, ]

    J1_i <- t(D_est_i) %*% V_est_i_inv %*% D_est_i
    J11 <- J11 + J1_i
    J2_i <- t(D_est_i) %*% V_est_i_inv

    Sigma1_i <- t(D_est_i) %*% AWA_est_i %*% (D_est_i)
    Sigma11 <- Sigma11 + Sigma1_i
    Sigma2_i <- t(D_est_i) %*% AWA_est_i

    VS_est_list <- c(VS_est_list, list(V_est_i_inv %*% S_est_i))
    AWA_est_list <- c(AWA_est_list, list(AWA_est_i))
    J2_list <- c(J2_list, list(J2_i))
    Sigma2_list <- c(Sigma2_list, list(Sigma2_i))
  }
  J11_inv <- chol2inv(chol(J11))
  JSigma11 <- J11_inv %*% Sigma11 %*% J11_inv
  res <- list(VS.est_list = VS_est_list,
              AWA.est_list = AWA_est_list,
              J2_list = J2_list,
              J11.inv = J11_inv,
              Sigma2_list = Sigma2_list,
              JSigma11 = JSigma11,
              mu.est = mu_est,
              V.est = V_est)
  return(res)
}