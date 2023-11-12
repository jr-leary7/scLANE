#' Truncated p-th power function (positive part).
#'
#' @name tp1
#' @param x A predictor variable value. Defaults to NULL.
#' @author Jakub Stoklosa and David I. Warton
#' @param t A specified knot value. Defaults to NULL.
#' @param p The \eqn{p^{th}} degree of the polynomial considered. Defaults to 1.
#' @return A vector of values that have been transformed using a truncated \eqn{p^{th}} power function (positive part) for a specified knot value.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @seealso \code{\link{tp2}}

tp1 <- function(x = NULL,
                t = NULL,
                p = 1) {
  res <- ((x - t)^p) * (x > t)
  return(res)
}
