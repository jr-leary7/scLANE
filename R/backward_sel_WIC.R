#' Backward selection function for MARGE - uses the Wald information criterion (WIC).
#'
#' @name backward_sel_WIC
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @author Jack Leary
#' @importFrom gamlss gamlss
#' @importFrom geeM geem
#' @importFrom MASS negative.binomial
#' @importFrom stats vcov
#' @param Y The response variable. Defaults to NULL.
#' @param B_new The model matrix. Defaults to NULL.
#' @param is.gee Is the model a GEE? Defaults to FALSE.
#' @param id.vec A vector of observation IDs that is necessary for fitting a GEE model. Defaults to NULL.
#' @param cor.structure The specified working correlation structure of the GEE model. Must be one of "independence", "ar1", or "exchangeable". Defaults to NULL.
#' @param sandwich.var Should the sandwich variance estimator be used instead of the model-based estimator? Default to FALSE.
#' @param theta.hat An initial estimate of \eqn{\hat{\theta}} used to fit the negative-binomial model when GEE mode is being used.
#' @return \code{backward_sel_WIC} returns the Wald statistic from the fitted model (the penalty is applied later on).
#' @references Stoklosa, J. Gibb, H. Warton, D.I. Fast forward selection for Generalized Estimating Equations With a Large Number of Predictor Variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.

backward_sel_WIC <- function(Y = NULL,
                             B_new = NULL,
                             is.gee = FALSE,
                             id.vec = NULL,
                             cor.structure = NULL,
                             theta.hat = NULL,
                             sandwich.var = FALSE) {
  # check inputs
  if (is.null(Y) || is.null(B_new)) { stop("Some inputs are missing from backward_sel_WIC().") }
  if (is.gee && is.null(id.vec)) { stop("GEEs require a vector of observation IDs in backward_sel_WIC().") }
  if (is.gee && is.null(cor.structure)) { stop("GEEs require a working correlation structure in backward_sel_WIC().") }
  if (is.gee && is.null(theta.hat)) { stop("In GEE mode you must provide an estimated value of theta to backward_sel_WIC().") }
  cor.structure <- tolower(cor.structure)
  if (is.gee) {
    fit <- geeM::geem(Y ~ B_new - 1,
                      id = id.vec,
                      corstr = cor.structure,
                      family = MASS::negative.binomial(theta.hat, link = "log"),
                      sandwich = sandwich.var)
    wald_stat <- unname(summary(fit)$wald.test[-1])^2
  } else {
    fit <- gamlss::gamlss(Y ~ B_new - 1,
                          family = "NBI",
                          trace = FALSE)
    vcov_mat <- try({ stats::vcov(fit, type = "all") }, silent = TRUE)
    if (inherits(vcov_mat, "try-error")) {
      covmat_unscaled <- chol2inv(fit$mu.qr$qr[1:(fit$mu.df - fit$mu.nl.df), 1:(fit$mu.df - fit$mu.nl.df), drop = FALSE])
      wald_stat <- unname(coef(fit) / sqrt(diag(covmat_unscaled))^2)[-1]
    } else {
      wald_stat <- unname((vcov_mat$coef / vcov_mat$se)[-c(1, length(vcov_mat$coef))]^2)
    }
  }
  return(wald_stat)
}
