#' Backward selection function for MARGE - uses the Wald information criterion (WIC).
#'
#' @name backward_sel_WIC
#' @param Y : the response variable. Defaults to NULL.
#' @param B_new : the model matrix. Defaults to NULL.
#' @param is.gee Is the model a GEE? Defaults to FALSE.
#' @param id.vec A vector of observation IDs that is necessary for fitting a GEE model. Defaults to NULL.
#' @param cor.structure The specified working correlation structure of the GEE model. Must be one of "independence", "ar1", or "exchangeable". Defaults to NULL.
#' @return \code{backward_sel_WIC} returns the Wald statistic from the fitted model (the penalty is applied later on).
#' @author Jakub Stoklosa and David I. Warton
#' @references Stoklosa, J. Gibb, H. Warton, D.I. Fast forward selection for Generalized Estimating Equations With a Large Number of Predictor Variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom gamlss gamlss
#' @importFrom geeM geem
#' @importFrom MASS negative.binomial
#' @seealso \code{\link{backward_sel}}

backward_sel_WIC <- function(Y = NULL, B_new = NULL, is.gee = FALSE, id.vec = NULL, cor.structure = NULL) {
  # check inputs
  if (is.null(Y) | is.null(B_new)) { stop("Some inputs are missing from backward_sel_WIC().") }
  if (is.gee & is.null(id.vec)) { stop("GEEs require a vector of observation IDs in backward_sel_WIC().") }
  if (is.gee & is.null(cor.structure)) { stop("GEEs require a working correlation structure in backward_sel_WIC().") }
  cor.structure <- tolower(cor.structure)
  if (is.gee) {
    fit <- geeM::geem(Y ~ B_new - 1,
                      id = id.vec,
                      corstr = cor.structure,
                      family = MASS::negative.binomial(1),
                      sandwich = TRUE)
    wald_stat <- (unname(summary(fit)$wald.test[-1]))^2
  } else {
    fit <- gamlss::gamlss(Y ~ B_new - 1, family = "NBI", trace = FALSE, )
    sink(tempfile())
    fit_sum <- summary(fit)
    sink()
    fit_sum_mat <- as.matrix(fit_sum)
    wald_stat <- ((fit_sum_mat[, 3])[-c(1, nrow(fit_sum_mat))])^2
  }
  return(wald_stat)
}
