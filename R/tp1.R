#' Truncated p-th power function (positive part).
#'
#' @name tp1
#' @param x : a vector of predictor variable values.
#' @param t : a specified knot value.
#' @param p : the pth degree of the polynomial considered.
#' @return \code{tp1} returns a vector of values that have been transformed using a truncated p-th power function (positive part) for a specified knot value.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @export
#' @seealso \code{\link{tp2}}
#' @examples
#' data(leptrine)
#' dat1 <- leptrine[[1]]
#' X_pred <- dat1[, 1]  # One predictor used.
#' tp1(X_pred, 1)  # Knot value set at x = 1.

tp1 <- function(x, t, p = 1) {
  res <- ((x - t)^p)*(x > t)
  return(res)
}
