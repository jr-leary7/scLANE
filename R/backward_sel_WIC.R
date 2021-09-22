#' Backward selection function for MARGE - uses the Wald information criterion (WIC).
#'
#' @name backward_sel_WIC
#' @param Y : the response variable.
#' @param B_new : the model matrix.
#' @return \code{backward_sel_WIC} returns the Wald statistic from the fitted model (the penalty is applied later on).
#' @author Jakub Stoklosa and David I. Warton
#' @references Stoklosa, J. Gibb, H. Warton, D.I. Fast forward selection for Generalized Estimating Equations With a Large Number of Predictor Variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @export
#' @seealso \code{\link{backward_sel}}

backward_sel_WIC <- function(Y, B_new, ...) {
  # check inputs
  if (any(sapply(c(Y, B_new), is.null))) stop("Some inputs are missing from backward_sel_WIC().")
  fit <- gamlss::gamlss(Y ~ B_new - 1, family = "NBI", trace = FALSE)
  sink(tempfile())
  fit_sum <- summary(fit)
  sink()
  fit_sum_mat <- as.matrix(fit_sum)
  wald_stat <- ((fit_sum_mat[, 3])[-c(1, nrow(fit_sum_mat))])^2
  return(wald_stat)
}
