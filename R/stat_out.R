#' Fits a linear regression model and calculates RSS/GCV measures (used for MARS linear models).
#'
#' @name stat_out
#' @param Y : the response variable.
#' @param B1 : the model matrix of predictor variables.
#' @param TSS : total sum of squares.
#' @param GCV.null : GCV value for the intercept model.
#' @param pen : the set/fixed penalty used for the GCV (the default is 2).
#' @details See the \code{earth} package for more details on the output measures calculated here.
#' @return \code{stat_out} returns a list of values, consisting of: RSS, RSSq1, GCV1 and GCVq1 values for the fitted model.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Milborrow, S. (2017a). Notes on the \code{earth} package. Package vignette. Available at: \url{http://127.0.0.1:31355/library/earth/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{earth}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{http://CRAN.R-project.org/package = earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats lm.fit fitted
#' @seealso \code{\link{stat_out_score_glm_null}}

stat_out <- function(Y = NULL,
                     B1 = NULL,
                     TSS = NULL,
                     GCV.null = NULL,
                     pen = 2) {
  # check inputs
  if (is.null(Y) || is.null(B1) || is.null(TSS) || is.null(GCV.null)) { stop("Some of the arguments to stat_out() are missing.") }
  if (GCV.null == 0) { stop("GCV.null argument to stat_out() cannot be 0.") }
  N <- length(Y)
  reg <- stats::lm.fit(B1, Y)
  if (any(is.na(reg$coef))) {
    RSS1 <- RSSq1 <- GCV1 <- GCVq1 <- NA  # return list of NAs if any \hat{beta} are NA
  } else {
    df1a <- ncol(B1) + pen * (ncol(B1) - 1) / 2  # This matches the earth() package, SAS and Friedman (1991) penalty.
    RSS1 <- sum((Y - stats::fitted(reg))^2)
    RSSq1 <- 1 - RSS1 / TSS
    GCV1 <- RSS1 / (N * (1 - (df1a) / N)^2)
    GCVq1 <- 1 - GCV1 / GCV.null
  }
  res <- list(RSS1 = RSS1,
              RSSq1 = round(RSSq1, 10),
              GCV1 = GCV1,
              GCVq1 = round(GCVq1, 10))
  return(res)
}
