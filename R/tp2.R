#' Truncated p-th power function (negative part).
#'
#' @name tp2
#' @author Jakub Stoklosa
#' @author David I. Warton
#' @param x A predictor variable value. Defaults to NULL.
#' @param t A specified knot value. Defaults to NULL.
#' @param p The \eqn{p^{th}} degree of the polynomial considered. Defaults to 1.
#' @return A vector of values that have been transformed using a truncated 
#' \eqn{p^{th}} power function (negative part) for a specified knot value.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. 
#' \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating 
#' equation approach to multivariate adaptive regression splines. 
#' \emph{Journal of Computational and Graphical Statistics}, 
#' \strong{27}, 245--253.
#' @seealso \code{\link{tp1}}

tp2 <- function(
        x = NULL,
        t = NULL,
        p = 1) {
    res <- ((t - x)^p) * (x < t)
    return(res)
}
