#' Backward selection function for MARS that uses the generalized cross validation criterion (GCV).
#'
#' @name backward_sel
#' @importFrom stats lm.fit fitted
#' @param Y : the response variable.
#' @param B_new : the model matrix.
#' @param pen : the set/fixed penalty used for the GCV. The default is 2.
#' @param GCV.null : GCV value for the intercept model. The default is 0.001.
#' @return \code{backward_sel} returns the GCV from the fitted model.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Milborrow, S. (2017a). Notes on the \code{earth} package. Package vignette. Available at: \url{http://127.0.0.1:31355/library/earth/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{earth}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{http://CRAN.R-project.org/package = earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats lm.fit fitted
#' @seealso \code{\link{backward_sel_WIC}}

backward_sel <- function(Y, B_new, pen = 2, GCV.null = 0.001) {
  # check inputs
  if (any(unlist(lapply(c(Y, B_new), is.null)))) stop("Some inputs are missing from backward_sel().")
  if (GCV.null == 0) stop("GCV.null in backward_sel() cannot be set to 0.")
  N <- length(Y)
  n_pred <- ncol(B_new) - 1
  GCV1 <- rep(NA, n_pred)
  for (j in seq(n_pred)) {
    B_new1 <- as.matrix(B_new[, -(j + 1)])
    n_pred1 <- ncol(B_new1)
    mod1 <- stats::lm.fit(B_new1, Y)
    RSS_back <- sum((Y - stats::fitted(mod1))^2)
    p <- n_pred1 + pen * (n_pred1 - 1) / 2  # This matches the earth() package, SAS and Friedman (1991) penalty.
    GCV1[j] <- 1 - (RSS_back / (N * (1 - p / N)^2)) / GCV.null
  }
  return(GCV1)
}
