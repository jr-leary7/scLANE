#' Truncated p-th power function (negative part).
#'
#' @name tp2
#' @param x : a predictor variable value.
#' @param t : a specified knot value.
#' @param p : the pth degree of the polynomial considered.
#' @return \code{tp2} returns a vector of values that have been transformed using a truncated p-th power function (negative part) for a specified knot value.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @seealso \code{\link{tp1}}
#' @examples
#' data(leptrine)
#' dat1 <- leptrine[[1]]
#' X_pred <- dat1[, 1]  # One predictor used.
#' tp2(X_pred, 1)  # Knot value set at x = 1.

tp2 <- function(x, t, p = 1) {
  res <- ((t - x)^p) * (x < t)
  return(res)
}
