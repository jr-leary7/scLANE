#' Truncates the predictor variable value to exclude extreme values in knots selection.
#'
#' @name max_span
#' @param X_pred : a vector of values for the predictor variable.
#' @param q : the number of predictors used.
#' @param alpha : see Friedman (1991) equation (45). The default is 0.05.
#' @details Note that this equation comes from Friedman (1991) equation (45).
#' @return \code{max_span} returns a vector of truncated predictor variable values.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @export
max_span <- function(X_pred = NULL, q = NULL, alpha = 0.05) {
  if (any(sapply(c(X_pred, q), is.null))) stop("Some inputs to max_span() are missing.")
  N <- length(unique(X_pred))
  x <- sort(unique(X_pred))
  maxspan <- round((3 - log2(alpha / q)))
  x_new <- x[-c(1:maxspan, floor(N - maxspan + 1):N)]
  if (length(x_new) > 0) {
    res <- x_new
  } else {
    res <- x
  }
  return(res)
}
